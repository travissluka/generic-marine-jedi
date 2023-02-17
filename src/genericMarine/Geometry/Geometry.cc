/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

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
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace genericMarine {

// ----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & conf,
                    const eckit::mpi::Comm & comm)
    : comm_(comm) {

  // create grid from configuration
  atlas::util::Config gridConfig(conf.getSubConfiguration("grid"));
  atlas::RegularLonLatGrid atlasRllGrid(gridConfig);
  functionSpace_ = atlas::functionspace::StructuredColumns(atlasRllGrid,
                    atlas::option::halo(0));
  functionSpaceIncHalo_ = atlas::functionspace::StructuredColumns(atlasRllGrid,
                          atlas::option::halo(0));
  // TODO(travis) put the halo back. Until then, this will only work on a single PE

  // load landmask
  if (conf.has("landmask.filename")) {
    loadLandMask(conf);
  }

  // calulate grid area
  // Temporary approximation solution, for a global
  // regular latlon grid, need to change if involved with other types of grid.
  const atlas::functionspace::StructuredColumns & fspace =
    static_cast<atlas::functionspace::StructuredColumns>(functionSpace());
  double dx = 2. * M_PI * atlas::util::DatumIFS::radius() / fspace.grid().nxmax();
  auto lonlat_data = atlas::array::make_view<double, 2>(functionSpace().lonlat());
  atlas::Field area = functionSpace().createField<double>(
      atlas::option::levels(1) | atlas::option::name("area"));
  auto area_data = atlas::array::make_view<double, 2>(area);
  for (int i=0; i < functionSpace().size(); i++) {
    area_data(i, 0) = dx*dx*cos(lonlat_data(i, 1)*M_PI/180.);
  }
  extraFields_.add(area);

  // vertical unit
  atlas::Field fld = fspace.createField<double>(
    atlas::option::levels(1) | atlas::option::name("vunit"));
  auto fld_data = atlas::array::make_view<double, 2>(fld);
  for (int i=0; i < functionSpace().size(); i++) {
    fld_data(i, 0) = 1.0;
  }
  extraFields_.add(fld);

  // add field for rossby radius
  if (conf.has("rossby radius file")) {
    readRossbyRadius(conf.getString("rossby radius file"));
  }
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
  auto fd = atlas::array::make_view<int, 2>(globalLandMask);

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
      for (int i = 0; i < lon; i++)
        fd(idx++, 0) = dataLandMask[j][i];
  }

  // scatter to the PEs
  atlas::Field fld = functionSpace().createField<int>(
                     atlas::option::levels(1) |
                     atlas::option::name("gmask"));
  extraFields_.add(fld);
  // TODO(travis) dangerous, don't do this?
  static_cast<atlas::functionspace::StructuredColumns>(functionSpace()).scatter(
    globalLandMask, extraFields_.field("gmask"));
}

// ----------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  const atlas::functionspace::StructuredColumns & fspace =
    static_cast<atlas::functionspace::StructuredColumns>(functionSpace());

  int ny = static_cast<int>(fspace.grid().ny());
  int nx = static_cast<int>(((atlas::RegularLonLatGrid&)(fspace.grid())).nx() );

  os << "Geometry: nx = " << nx << ", ny = " << ny << std::endl;

  size_t nMaskedLand = 0;
  size_t nUnmaskedOcean = 0;
  const int nSize = functionSpace().size();
  auto fd = atlas::array::make_view<int, 2>(extraFields_.field("gmask"));
  for (size_t j = 0; j < nSize; j++) {
    if (fd(j, 0) == 1) {
      nUnmaskedOcean += 1;
    } else if (fd(j, 0) == 0) {
      nMaskedLand += 1;
    } else {
      util::abor1_cpp("Geometry::print(), landmask neither 1 nor 0.",
        __FILE__, __LINE__);
    }
  }
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
  for (int i=0; i < functionSpace().size(); i++) {
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

  return dstField;
}

// ----------------------------------------------------------------------------
void Geometry::latlon(std::vector<double> &lats, std::vector<double> & lons,
                      const bool halo) const {
  const atlas::functionspace::StructuredColumns* fspace;
  if (halo) {
    fspace = &functionSpaceIncHalo_;
  } else {
    fspace = &functionSpace_;
  }

  auto lonlat = atlas::array::make_view<double, 2>(fspace->lonlat());
  auto ngrid = lonlat.shape<0>();
  for (size_t i=0; i < ngrid; i++) {
    // TODO(travis) don't do this, how should the halos work at the poles?
    auto lat = lonlat(i, 1);
    if (lat > 90) lat = 90;
    if (lat < -90) lat = -90;
    lats.push_back(lat);
    lons.push_back(lonlat(i, 0));
  }
}

// ----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> lvls(vars.size(), 1);
  // TODO(travis) get the actual number of levels
  return lvls;
}

// ----------------------------------------------------------------------------
}  // namespace genericMarine
