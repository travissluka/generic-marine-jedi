/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <fstream>
#include <vector>
#include "netcdf"

#include "genericMarine/Geometry/Geometry.h"

#include "eckit/container/KDTree.h"
#include "eckit/config/Configuration.h"

#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/option.h"
#include "atlas/functionspace.h"
#include "atlas/util/Config.h"

#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace genericMarine {

// ----------------------------------------------------------------------------

const int GEOM_HALO_SIZE = 1;

// ----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & conf, const eckit::mpi::Comm & comm)
    : comm_(comm), extraFields_() {
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // create grid from configuration
  // NOTE: we have to use "checkerboard" instead of the default so that the
  // poles aren't placed completely on a single PE, which breaks the interpolation
  // at the moment due to a bug with redundant points.
  atlas::util::Config gridConfig(conf.getSubConfiguration("grid"));
  atlas::RegularLonLatGrid atlasRllGrid(gridConfig);
  functionSpace_ = atlas::functionspace::StructuredColumns(atlasRllGrid,
                      atlas::grid::Partitioner("checkerboard"),
                      atlas::option::halo(GEOM_HALO_SIZE) );
  atlas::functionspace::StructuredColumns fs(functionSpace_);
  auto vGhost = atlas::array::make_view<int, 1>(fs.ghost());
  auto vLonlat = atlas::array::make_view<double, 2>(fs.lonlat());

  // debugging information about the grid
  #define LL_INIT {9e9, -9e9}
  double lat[2] = LL_INIT, lon[2] = LL_INIT, hLat[2] = LL_INIT, hLon[2] = LL_INIT;
  for (int i = 0; i < functionSpace_.size(); i++) {
    hLon[0] = std::min(hLon[0], vLonlat(i, 0)); hLon[1] = std::max(hLon[1], vLonlat(i, 0));
    hLat[0] = std::min(hLat[0], vLonlat(i, 1)); hLat[1] = std::max(hLat[1], vLonlat(i, 1));
    if (vGhost(i)) continue;

    lon[0] = std::min(lon[0], vLonlat(i, 0)); lon[1] = std::max(lon[1], vLonlat(i, 0));
    lat[0] = std::min(lat[0], vLonlat(i, 1)); lat[1] = std::max(lat[1], vLonlat(i, 1));
  }
  oops::Log::debug() << "grid      (Lat)/(Lon): (" << lat[0] << ", " << lat[1] << ") / ("
                     << lon[0] << " , " << lon[1] << ")"<< std::endl;
  oops::Log::debug() << "grid halo (Lat)/(Lon): (" << hLat[0] << ", " << hLat[1] << ") / ("
                     << hLon[0] << " , " << hLon[1] << ")"<< std::endl;

  // load landmask, after this call the fields "mask" (floating point) and "gmask" (integer)
  // will be added. We need both because bump and the interpolation have different requirements
  if (conf.has("landmask.filename")) {
    loadLandMask(conf);
  }

  // calulate grid area
  // Temporary approximation solution, for a global
  // regular latlon grid, need to change if involved with other types of grid.
  double dx = 2. * M_PI * atlas::util::DatumIFS::radius() / fs.grid().nxmax();
  atlas::Field area = functionSpace().createField<double>(
      atlas::option::levels(1) | atlas::option::name("area"));
  auto vArea = atlas::array::make_view<double, 2>(area);
  for (int i=0; i < functionSpace().size(); i++) {
    double lat = vLonlat(i, 1);
    if (lat > 90) lat = 180 - lat;
    if (lat < -90) lat = -180 - lat;
    vArea(i, 0) = dx*dx*cos(lat*M_PI/180.);
  }
  extraFields_.add(area);

  // vertical unit
  atlas::Field vunit = fs.createField<double>(
    atlas::option::levels(1) | atlas::option::name("vunit"));
  auto vVunit = atlas::array::make_view<double, 2>(vunit);
  vVunit.assign(1.0);
  extraFields_.add(vunit);

  // add field for rossby radius
  if (conf.has("rossby radius file")) {
    readRossbyRadius(conf.getString("rossby radius file"));
  }

  // halo mask (needed for SABER)
  // NOTE: this has to be done AFTER the halo exchange, otherwise it would be all 1 !
  atlas::Field hmask = functionSpace().createField<int>(
      atlas::option::levels(1) | atlas::option::name("hmask"));
  auto vHmask = atlas::array::make_view<int, 2>(hmask);
  vHmask.assign(0);
  for (int i=0; i < fs.size(); i++) {
    if (vGhost(i)) continue;
    vHmask(i, 0) = 1;
  }
  extraFields_.add(hmask);
}

// ----------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_) {
  ASSERT(1 == 2);
}

// ----------------------------------------------------------------------------

Geometry::~Geometry() {}

// ----------------------------------------------------------------------------

void Geometry::loadLandMask(const eckit::Configuration &conf) {
  // use an globalLandMask to read the data on root PE only.
  atlas::Field globalLandMask = functionSpace().createField<int>(
                                atlas::option::levels(1) |
                                atlas::option::name("gmask") |
                                atlas::option::global());

  const int size = functionSpace().size();
  auto vGlobalLandMask = atlas::array::make_view<int, 2>(globalLandMask);

  // read file only on the root PE.
  if (globalLandMask.size() != 0) {
    int lat = 0, lon = 0;
    std::string filename;

    if (!conf.get("landmask.filename", filename))
      util::abor1_cpp("Geometry::loadLandMask(), Get filename failed.",
        __FILE__, __LINE__);
    else
      oops::Log::info() << "In Geometry::loadLandMask(), filename = "
                        << filename << std::endl;

    // Open netCDF file
    netCDF::NcFile file(filename.c_str(), netCDF::NcFile::read);
    if (file.isNull())
      util::abor1_cpp("Geometry::loadLandMask(), Create netCDF file failed.",
        __FILE__, __LINE__);

    // get file dimensions
    lat = static_cast<int>(file.getDim("lat").getSize());
    lon = static_cast<int>(file.getDim("lon").getSize());

    // get landmask data
    netCDF::NcVar varLandMask;
    varLandMask = file.getVar("landmask");
    if (varLandMask.isNull())
      util::abor1_cpp("Get var landmask failed.", __FILE__, __LINE__);

    int dataLandMask[lat][lon];
    varLandMask.getVar(dataLandMask);

    // TODO(someone) the netcdf lat dimension is likely inverted compared to
    // the  atlas grid. This should be explicitly checked.
    int idx = 0;
    for (int j = lat-1; j >= 0; j--)
      for (int i = 0; i < lon; i++) {
        vGlobalLandMask(idx++, 0) = dataLandMask[j][i];
      }
  }

  // scatter to the PEs
  // gmask is needed for SABER
  atlas::Field gmask = functionSpace().createField<int>(
                     atlas::option::levels(1) |
                     atlas::option::name("gmask"));
  functionSpace_.scatter(globalLandMask, gmask);
  functionSpace().haloExchange(gmask);
  extraFields_.add(gmask);

  // create a floating point version
  // mask is needed for OOPS interpolation
  atlas::Field mask = functionSpace().createField<double>(
                        atlas::option::levels(1) |
                        atlas::option::name("mask"));
  auto vGmask = atlas::array::make_view<int, 2>(gmask);
  auto vMask = atlas::array::make_view<double, 2>(mask);
  for (int j = 0; j < size; j++)
    vMask(j, 0) = static_cast<double>(vGmask(j, 0));
  extraFields_.add(mask);
}

// ----------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  // get some variables we'll need later
  atlas::functionspace::StructuredColumns fspace(functionSpace_);
  atlas::StructuredGrid grid = fspace.grid();
  atlas::idx_t ny = grid.ny();
  atlas::idx_t nx = grid.nxmax();

  // global size
  os << "Geometry: nx = " << nx << ", ny = " << ny << std::endl;

  // count grid points on local PE
  size_t nMaskedLand = 0;
  size_t nUnmaskedOcean = 0;
  const int nSize = functionSpace().size();
  auto fd = atlas::array::make_view<int, 2>(extraFields_.field("gmask"));
  auto ghost = atlas::array::make_view<int, 1>(functionSpace().ghost());
  for (size_t j = 0; j < nSize; j++) {
    // don't count halo (ghost) points
    if (ghost(j) != 0) continue;

    if (fd(j, 0) == 1) {
      nUnmaskedOcean += 1;
    } else if (fd(j, 0) == 0) {
      nMaskedLand += 1;
    } else {
      util::abor1_cpp("Geometry::print(), landmask neither 1 nor 0.",
        __FILE__, __LINE__);
    }
  }

  // gather results from all PEs
  oops::mpi::world().allReduceInPlace(nMaskedLand, eckit::mpi::Operation::SUM);
  oops::mpi::world().allReduceInPlace(nUnmaskedOcean, eckit::mpi::Operation::SUM);

  os << "Geometry: # of unmasked ocean grid = " << nUnmaskedOcean
     << ", # of masked land grid = " << nMaskedLand << std::endl;
}

// ----------------------------------------------------------------------------

void Geometry::readRossbyRadius(const std::string & filename) {
  std::ifstream infile(filename);
  std::vector<eckit::geometry::Point2> lonlat;
  std::vector<double> vals;
  double lat, lon, x, val;

  while (infile >> lat >> lon >> x >> val) {
    lonlat.push_back(eckit::geometry::Point2(lon, lat));
    vals.push_back(val*1.0e3);
  }

  atlas::Field field = interpToGeom(lonlat, vals);
  field.rename("rossby_radius");
  extraFields_.add(field);
}

// ----------------------------------------------------------------------------

atlas::Field Geometry::interpToGeom(
  const std::vector<eckit::geometry::Point2> & srcLonLat,
  const std::vector<double> & srcVal) const
{
  // TODO(travis) replace this with atlas interpolation
  // Interpolate the values from the given lat/lons onto the grid that is
  // represented by this geometry. Note that this assumes each PE is
  // presenting an identical copy of srcLonLat and srcVal.
  struct TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef double                  Payload;
  };
  typedef eckit::KDTreeMemory<TreeTrait> KDTree;
  const int maxSearchPoints = 4;

  // Create a KD tree for fast lookup
  std::vector<typename KDTree::Value> srcPoints;
  for (int i = 0; i < srcVal.size(); i++) {
    KDTree::PointType xyz;
    atlas::util::Earth::convertSphericalToCartesian(srcLonLat[i], xyz);
    srcPoints.push_back(KDTree::Value(xyz, srcVal[i]) );
  }
  KDTree kd;
  kd.build(srcPoints.begin(), srcPoints.end());

  // Interpolate (inverse distance weighted)
  atlas::Field dstField = functionSpace().createField<double>(
    atlas::option::levels(1));
  auto dstView = atlas::array::make_view<double, 2>(dstField);
  auto dstLonLat = atlas::array::make_view<double, 2>(functionSpace().lonlat());
  auto vGhost = atlas::array::make_view<int, 1>(functionSpace_.ghost());
  for (int i=0; i < functionSpace().size(); i++) {
    if (vGhost(i)) continue;

    eckit::geometry::Point2 dstPoint({dstLonLat(i, 0), dstLonLat(i, 1)});
    eckit::geometry::Point3 dstPoint3D;
    atlas::util::Earth::convertSphericalToCartesian(dstPoint, dstPoint3D);
    auto points = kd.kNearestNeighbours(dstPoint3D, maxSearchPoints);
    double sumDist = 0.0;
    double sumDistVal = 0.0;
    for ( int n = 0; n < points.size(); n++ ) {
      if ( points[n].distance() < 1.0e-6 ) {
        sumDist = 1.0;
        sumDistVal = points[n].payload();
        break;
      }
      double w = 1.0 / (points[n].distance()*points[n].distance());
      sumDist += w;
      sumDistVal += w*points[n].payload();
    }

    dstView(i, 0) = sumDistVal / sumDist;
  }

  functionSpace_.haloExchange(dstField);
  return dstField;
}

// ----------------------------------------------------------------------------
void Geometry::latlon(std::vector<double> &lats, std::vector<double> & lons,
                      const bool halo) const {
  atlas::functionspace::StructuredColumns fs(functionSpace_);

  size_t size = halo ? fs.sizeHalo() : fs.sizeOwned();
  lats.resize(size);
  lons.resize(size);

  // get the list of lat/lons
  auto vLonlat = atlas::array::make_view<double, 2>(fs.lonlat());
  auto vGhost = atlas::array::make_view<int, 1>(fs.ghost());
  int j = 0;
  for (size_t i=0; i < fs.size(); i++) {
    // skip if a halo point, and if we don't want halo points
    if (!halo && vGhost(i)) continue;

    double lon = vLonlat(i, 0);
    double lat = vLonlat(i, 1);

    // TODO(travis) don't do this, how should the halos work at the poles?
    // this is a hack, and is wrong. Wrap the halo points beyond the poles.
    if (vGhost(i) && lat < -90 ) { lat = -180. - lat; lon += 180.; }
    if (vGhost(i) && lat > 90 ) { lat = 180. - lat; lon += 180.; }
    // if (lon > 180) lon -= 360;

    lats[j] = lat;
    lons[j++] = lon;
  }
  ASSERT(j == lats.size());
}

// ----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> lvls(vars.size(), 1);
  // TODO(travis) get the actual number of levels
  return lvls;
}

// ----------------------------------------------------------------------------
}  // namespace genericMarine
