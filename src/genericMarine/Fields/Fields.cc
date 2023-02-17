/*
 * (C) Copyright 2020-2023 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <limits>
#include <string>

#include "netcdf"

#include "genericMarine/Fields/Fields.h"
#include "genericMarine/Geometry/Geometry.h"
#include "genericMarine/State/State.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/option.h"

#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

using atlas::array::make_view;

namespace genericMarine {

// ----------------------------------------------------------------------------

Fields::Fields(const Geometry & geom, const oops::Variables & vars,
                const util::DateTime & vt)
  : atlasFieldSet_(), geom_(geom), missing_(util::missingValue(this->missing_)),
    time_(vt), vars_(vars) {
  // the constructor that gets called by everything (all State and Increment
  //  constructors ultimately end up here)
  updateFields(vars);
}

// ----------------------------------------------------------------------------

Fields::Fields(const Fields & other)
  : Fields(other.geom_, other.vars_, other.time_) {
  // copy data from object other
  *this = other;
}

// ----------------------------------------------------------------------------

void Fields::updateFields(const oops::Variables & vars) {
  atlas::FieldSet fset;
  for (int v = 0; v < vars.size(); v++) {
    if (atlasFieldSet_.has(vars[v])) {
      // field already exists, copy over
      fset.add(atlasFieldSet_.field(vars[v]));
    } else {
      // field does not exist, create
      atlas::Field fld = geom_.functionSpace().createField<double>(
                          atlas::option::levels(1) |
                          atlas::option::name(vars[v]));
      auto fd = make_view<double, 2>(fld);
      fd.assign(0.0);
      fset.add(fld);
    }
  }
  atlasFieldSet_ = fset;
  vars_ = vars;
}

// ----------------------------------------------------------------------------

Fields & Fields::operator =(const Fields & other) {
  time_ = other.time_;

  updateFields(other.vars_);

  const int size = geom_.functionSpace().size();
  for (int v = 0; v < vars_.size(); v++) {
    std::string name = vars_[v];
    auto fd       = make_view<double, 2>(atlasFieldSet_.field(name));
    auto fd_other = make_view<double, 2>(other.atlasFieldSet_.field(name));

    for (int j = 0; j < size; j++)
      fd(j, 0) = fd_other(j, 0);
  }
  return *this;
}

// ----------------------------------------------------------------------------

Fields & Fields::operator+=(const Fields &other) {
  const int size = geom_.functionSpace().size();
  for (int v = 0; v < vars_.size(); v++) {
    std::string name = vars_[v];
    ASSERT(other.atlasFieldSet_.has(name));
    auto fd       = make_view<double, 2>(atlasFieldSet_.field(name));
    auto fd_other = make_view<double, 2>(other.atlasFieldSet_.field(name));

    for (int j = 0; j < size; j++) {
      if (fd(j, 0) == missing_ || fd_other(j, 0) == other.missing_)
        fd(j, 0) = missing_;
      else
        fd(j, 0) += fd_other(j, 0);
    }
  }
  return *this;
}

// ----------------------------------------------------------------------------

void Fields::accumul(const double &zz, const Fields &rhs) {
  const size_t size = geom_.functionSpace().size();

  for (int v = 0; v < vars_.size(); v++) {
    std::string name = vars_[v];
    ASSERT(rhs.atlasFieldSet_.has(name));
    auto fd = make_view<double, 2>(atlasFieldSet_.field(name));
    auto fd_rhs = make_view<double, 2>(rhs.atlasFieldSet_.field(name));

    for (size_t i = 0; i < size; i++) {
      if (fd(i, 0) == missing_ || fd_rhs(i, 0) == missing_)
        fd(i, 0) = missing_;
      else
        fd(i, 0) +=  zz*fd_rhs(i, 0);
    }
  }
}

// ----------------------------------------------------------------------------

double Fields::norm() const {
  const int size = geom_.functionSpace().size();
  int nValid = 0;
  double norm = 0.0, s = 0.0;

  for (int v = 0; v < vars_.size(); v++) {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(v));

    for (int i = 0; i < size; i++) {
      if (fd(i, 0) != missing_) {
        nValid += 1;
        s += fd(i, 0)*fd(i, 0);
      }
    }
  }

  // sum results across PEs
  oops::mpi::world().allReduceInPlace(nValid, eckit::mpi::Operation::SUM);
  oops::mpi::world().allReduceInPlace(s, eckit::mpi::Operation::SUM);

  if (nValid == 0)
    norm = 0.0;
  else
    norm = sqrt(s/(1.0*nValid));

  return norm;
}

// ----------------------------------------------------------------------------

void Fields::zero() {
  const int size = geom_.functionSpace().size();
  for (int v = 0; v < vars_.size(); v++) {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(v));
    fd.assign(0.0);
  }
}

// ----------------------------------------------------------------------------

void Fields::read(const eckit::Configuration & conf) {
    // create a global field valid on the root PE
  atlas::Field globalSst = geom_.functionSpace().createField<double>(
                       atlas::option::levels(1) |
                       atlas::option::global());

  // following code block should execute on the root PE only
  if ( globalSst.size() != 0 ) {
    int time = 0, lon = 0, lat = 0;
    std::string filename;

    auto fd = make_view<double, 2>(globalSst);

    const atlas::functionspace::StructuredColumns & fspace =
      static_cast<atlas::functionspace::StructuredColumns>(geom_.functionSpace());

    // get filename
    if (!conf.get("filename", filename))
      util::abor1_cpp("Fields::read(), Get filename failed.",
        __FILE__, __LINE__);

    // open netCDF file
    netCDF::NcFile file(filename.c_str(), netCDF::NcFile::read);
    if (file.isNull())
      util::abor1_cpp("Fields::read(), Create netCDF file failed.",
        __FILE__, __LINE__);

    // get file dimensions
    time = static_cast<int>(file.getDim("time").getSize());
    lon  = static_cast<int>(file.getDim("lon").getSize());
    lat  = static_cast<int>(file.getDim("lat").getSize());
    {
      // sanity check
      int ny = static_cast<int>(fspace.grid().ny());
      int nx = static_cast<int>(((atlas::RegularLonLatGrid&)(fspace.grid())).nx() );
      if (time != 1 || lat != ny || lon != nx) {
        util::abor1_cpp("Fields::read(), lat!=ny or lon!=nx",
        __FILE__, __LINE__);
      }
    }

    // get sst data
    netCDF::NcVar sstVar;
    sstVar = file.getVar("sst");
    if (sstVar.isNull())
      util::abor1_cpp("Get sst var failed.", __FILE__, __LINE__);
    float  sstData[lat][lon];
    sstVar.getVar(sstData);  // if used double, read-in data would be wrong.

    // mask missing values
    const double epsilon = 1.0e-6;
    const double missing_nc = -32768.0;
    bool isKelvin = conf.getBool("kelvin", false);
    for (int j = 0; j < lat; j++)
      for (int i = 0; i < lon; i++)
        if (abs(sstData[j][i]-(missing_nc)) < epsilon) {
          sstData[j][i] = missing_;
          // TODO(someone) missing values that aren't a part of the landmask
          // should be filled in instead
        } else {
          // Kelvin to Celsius which JEDI use internally, will check if the
          // units is Kelvin or Celsius in the future
          if (isKelvin)
            sstData[j][i] -= 273.15;
        }

    // float to double
    int idx = 0;
    for (int j = lat-1; j >= 0; j--)
      for (int i = 0; i < lon; i++)
        fd(idx++, 0) = static_cast<double>(sstData[j][i]);
  }

  // scatter to the PEs
  // TODO(travis) dangerous, don't do this?
  static_cast<atlas::functionspace::StructuredColumns>(geom_.functionSpace()).scatter(
    globalSst, atlasFieldSet_.field("sea_surface_temperature"));

  // apply mask from read in landmask
  if ( geom_.extraFields().has("gmask") ) {
    atlas::Field mask_field = geom_.extraFields()["gmask"];
    auto mask = make_view<int, 2>(mask_field);
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    for (int i = 0; i < mask.size(); i++) {
      if (mask(i, 0) == 0)
        fd(i, 0) = missing_;
    }
  }

  // done, do halo exchange
  atlasFieldSet_.haloExchange();
}

// ----------------------------------------------------------------------------

void Fields::write(const eckit::Configuration & conf) const {
  const atlas::functionspace::StructuredColumns & fspace =
    static_cast<atlas::functionspace::StructuredColumns>(geom_.functionSpace());

  // gather from the PEs
  atlas::Field globalSst = fspace.createField<double>(
    atlas::option::levels(1) |
    atlas::option::global());
  fspace.gather(atlasFieldSet_.field("sea_surface_temperature"), globalSst);

  // The following code block should execute on the root PE only
  if ( globalSst.size() != 0 ) {
    int lat, lon, time = 1;
    std::string filename;

     oops::Log::info() << "In Fields::write(), conf = " << conf << std::endl;

    // get filename
    if (!conf.get("filename", filename)) {
        util::abor1_cpp("Fields::write(), Get filename failed.",
                        __FILE__, __LINE__);
    } else {
      oops::Log::info() << "Fields::write(), filename=" << filename
                        << std::endl;
    }

    // create netCDF file
    netCDF::NcFile file(filename.c_str(), netCDF::NcFile::replace);
    if (file.isNull())
      util::abor1_cpp("Fields::write(), Create netCDF file failed.",
                      __FILE__, __LINE__);

    // define dims
    lat = fspace.grid().ny();
    lon = ((atlas::RegularLonLatGrid)(fspace.grid())).nx();

    // unlimited dim if without size parameter, then it'll be 0,
    // what about the size?
    netCDF::NcDim timeDim = file.addDim("time", 1);
    netCDF::NcDim latDim  = file.addDim("lat" , lat);
    netCDF::NcDim lonDim  = file.addDim("lon" , lon);
    if (timeDim.isNull() || latDim.isNull() || lonDim.isNull())
      util::abor1_cpp("Fields::write(), Define dims failed.",
                      __FILE__, __LINE__);

    std::vector<netCDF::NcDim> dims;
    dims.push_back(timeDim);
    dims.push_back(latDim);
    dims.push_back(lonDim);

    // inside use double, read/write use float
    // to make it consistent with the netCDF files.
    netCDF::NcVar sstVar = file.addVar(std::string("sst"),
                                      netCDF::ncFloat, dims);

    // define units atts for data vars
    const float fillvalue = -32768.0;
    sstVar.putAtt("units", "K");
    sstVar.putAtt("_FillValue", netCDF::NcFloat(), fillvalue);
    sstVar.putAtt("missing_value", netCDF::NcFloat(), fillvalue);

    // write data to the file
    auto fd = atlas::array::make_view<double, 2>(globalSst);
    float sstData[time][lat][lon];
    bool isKelvin = conf.getBool("kelvin", false);
    int idx = 0;
    for (int j = lat-1; j >= 0; j--)
      for (int i = 0; i < lon; i++) {
        if (fd(idx, 0) == missing_) {
          sstData[0][j][i] = fillvalue;
        } else {
          // doulbe to float, also convert JEDI Celsius to Kelvin, in the
          // future it should be able to handle both Kelvin and Celsius.
          if (isKelvin)
            sstData[0][j][i] = static_cast<float>(fd(idx, 0)) + 273.15;
          else
            sstData[0][j][i] = static_cast<float>(fd(idx, 0));
        }
        idx++;
      }

    sstVar.putVar(sstData);

    oops::Log::info() << "Fields::write(), Successfully write data to file!"
                      << std::endl;
  }
}

// ----------------------------------------------------------------------------

void Fields::toFieldSet(atlas::FieldSet & fset) const {
  const int size = geom_.functionSpace().size();
  for (int v = 0; v < vars_.size(); v++) {
    std::string name = vars_[v];
    ASSERT(atlasFieldSet_.has(name));

    atlas::Field fld = geom_.functionSpace().createField<double>(
                        atlas::option::levels(1) |
                        atlas::option::name(name));
    fld.metadata().set("interp_type", "default");

    auto fd  = make_view<double, 2>(fld);
    auto fd2 = make_view<double, 2>(atlasFieldSet_.field(name));
    for (int j = 0; j < size; j++)
      fd(j, 0) = fd2(j, 0);

    fset.add(fld);
  }
}

// ----------------------------------------------------------------------------

void Fields::toFieldSetAD(const atlas::FieldSet & fset) {
  const int size = geom_.functionSpace().size();
  for (int v = 0; v < vars_.size(); v++) {
    std::string name = vars_[v];
    ASSERT(fset.has(name));

    auto fd    = make_view<double, 2>(atlasFieldSet_.field(name));
    auto fd_in = make_view<double, 2>(fset.field(name));

    for (int j = 0; j < size; j++) {
      if (fd(j, 0) == missing_ || fd_in(j, 0) == missing_)
        fd(j, 0) = missing_;
      else
        fd(j, 0) += fd_in(j, 0);
    }
  }

  // TODO(travis) adjoint halo exchange
  // atlasFieldSet_.adjointHaloExchange();
}

// ----------------------------------------------------------------------------

void Fields::fromFieldSet(const atlas::FieldSet & fset) {
  const int size = geom_.functionSpace().size();
  for (int v = 0; v < vars_.size(); v++) {
    std::string name = vars_[v];
    ASSERT(fset.has(name));

    auto fd    = make_view<double, 2>(atlasFieldSet_.field(name));
    auto fd_in = make_view<double, 2>(fset.field(name));

    for (int j = 0; j < size; j++) {
      fd(j, 0) = fd_in(j, 0);
    }
  }
}

// ----------------------------------------------------------------------------

void Fields::print(std::ostream & os) const {
  const int size = geom_.functionSpace().size();
  for (int v = 0; v < vars_.size(); v++) {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(v));
    double mean = 0.0, sum = 0.0,
          min = std::numeric_limits<double>::max(),
          max = std::numeric_limits<double>::min();
    int nValid = 0;

    for (int i = 0; i < size; i++)
      if (fd(i, 0) != missing_) {
        if (fd(i, 0) < min) min = fd(i, 0);
        if (fd(i, 0) > max) max = fd(i, 0);

        sum += fd(i, 0);
        nValid++;
      }

    // gather results across PEs
    oops::mpi::world().allReduceInPlace(nValid, eckit::mpi::Operation::SUM);
    oops::mpi::world().allReduceInPlace(sum, eckit::mpi::Operation::SUM);
    oops::mpi::world().allReduceInPlace(min, eckit::mpi::Operation::MIN);
    oops::mpi::world().allReduceInPlace(max, eckit::mpi::Operation::MAX);

    if (nValid == 0) {
      mean = 0.0;
      oops::Log::debug() << "Field::print(), nValid == 0!" << std::endl;
    } else {
      mean = sum / (1.0*nValid);
    }

    os << "min = " << min << ", max = " << max << ", mean = " << mean
      << std::endl;
  }
}
// ----------------------------------------------------------------------------

}  // namespace genericMarine
