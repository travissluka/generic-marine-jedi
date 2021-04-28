/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "netcdf"

#include "umdsst/Geometry/Geometry.h"

#include "eckit/config/Configuration.h"

#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

using atlas::array::make_view;

namespace umdsst {

// ----------------------------------------------------------------------------

  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm) : comm_(comm) {
    atlas::util::Config gridConfig(conf.getSubConfiguration("grid"));
    atlas::RegularLonLatGrid atlasRllGrid(gridConfig);

    atlasFunctionSpace_.reset(
      new atlas::functionspace::StructuredColumns(atlasRllGrid,
      atlas::option::halo(0)));

    atlasFieldSet_.reset(new atlas::FieldSet());
    atlasFieldSet_->add(atlasFunctionSpace_->lonlat());

    if (conf.has("landmask.filename")) {
      oops::Log::debug() << "Geometry::Geometry(), before loading landmask."
                        << std::endl;
      loadLandMask(conf);

      oops::Log::debug() << "Geometry::Geometry(), after loading landmask."
                        << std::endl;
    }

    // add Field area
    // Ligang: From Travis, Temporary approximation solution, for a global
    // regular latlon grid, need to change if involved with other types of grid.
    double dx = 2. * M_PI * atlas::util::DatumIFS::radius()
        / atlasFunctionSpace()->grid().nxmax();
    auto lonlat_data = make_view<double, 2>(atlasFunctionSpace_->lonlat());
    atlas::Field area = atlasFunctionSpace_->createField<double>(
        atlas::option::levels(1) | atlas::option::name("area"));
    auto area_data = make_view<double, 2>(area);
    for (int i=0; i < atlasFunctionSpace_->size(); i++) {
      area_data(i, 0) = dx*dx*cos(lonlat_data(i, 1)*M_PI/180.);
    }
    atlasFieldSet_->add(area);
  }

// ----------------------------------------------------------------------------

  Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
    atlasFunctionSpace_.reset(new
      atlas::functionspace::StructuredColumns(other.atlasFunctionSpace_->grid(),
      atlas::option::halo(0)));

    atlasFieldSet_.reset(new atlas::FieldSet());
//  atlasFieldSet_->add(atlasFunctionSpace_->lonlat());
    for (int i = 0; i < other.atlasFieldSet()->size(); i++) {
      atlasFieldSet_->add((*other.atlasFieldSet())[i]);
    }
  }

// ----------------------------------------------------------------------------

  Geometry::~Geometry() {}

// ----------------------------------------------------------------------------

  void Geometry::loadLandMask(const eckit::Configuration &conf) {
    // use an globalLandMask to read the data on root PE only.
    atlas::Field globalLandMask = atlasFunctionSpace_->createField<int>(
                                  atlas::option::levels(1) |
                                  atlas::option::name("gmask") |
                                  atlas::option::global());
    auto fd = make_view<int, 2>(globalLandMask);

    // Ligang: read file only on the root PE.
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

    atlas::Field fld = atlasFunctionSpace_->createField<int>(
                       atlas::option::levels(1) |
                       atlas::option::name("gmask"));
    atlasFieldSet_->add(fld);

    atlasFunctionSpace_->scatter(globalLandMask,
                                 atlasFieldSet_->field("gmask"));
  }

// ----------------------------------------------------------------------------

  void Geometry::print(std::ostream & os) const {
    int nx, ny;
    ny = static_cast<int>(atlasFunctionSpace()->grid().ny());
    nx = static_cast<int>(
      ((atlas::RegularLonLatGrid&)(atlasFunctionSpace()->grid())).nx() );

    os << "Geometry: nx = " << nx << ", ny = " << ny << std::endl;

    int nMaskedLand = 0, nUnmaskedOcean = 0,
        nSize = atlasFunctionSpace_->size();
    auto fd = make_view<int, 2>(atlasFieldSet_->field("gmask"));
    for (int j = 0; j < nSize; j++) {
      if (fd(j, 0) == 1)
        nUnmaskedOcean += 1;
      else if (fd(j, 0) == 0)
        nMaskedLand += 1;
      else
        util::abor1_cpp("Geometry::print(), landmask neither 1 nor 0.",
          __FILE__, __LINE__);
    }
    os << "Geometry: # of unmasked ocean grid = " << nUnmaskedOcean
       << ", # of masked land grid = " << nMaskedLand << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace umdsst
