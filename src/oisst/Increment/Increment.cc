/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>
#include <limits>

#include "netcdf"

#include "oisst/Geometry/Geometry.h"
#include "oisst/Increment/Increment.h"
#include "oisst/State/State.h"

#include "eckit/config/Configuration.h"

#include "atlas/field.h"
#include "atlas/array.h"
#include "atlas/option.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

// using namespace netCDF;
using atlas::array::make_view;
using atlas::option::name;


namespace oisst {

// ----------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom,
                       const oops::Variables & vars,
                       const util::DateTime & vt)
    : geom_(new Geometry(geom)), time_(vt), vars_(vars) {

    if (vars_.size() != 1) {
      util::abor1_cpp("Increment::Increment(), vars_.size() != 1.",
                      __FILE__, __LINE__);
    }

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<double>(
                         name(var));

      auto fd = make_view<double, 1>(fld);
      fd.assign(0.0);

      atlasFieldSet_->add(fld);
    }
  }

// ----------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom, const Increment & other)
    : geom_(new Geometry(geom)), time_(other.time_), vars_(other.vars_) {
    // it will normally be used for interpolation and change resolution.

    if (vars_.size() != 1) {
      util::abor1_cpp("Increment::Increment(), vars_.size() != 1.",
                       __FILE__, __LINE__);
    }

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<double>(
                         name(var));
      atlasFieldSet_->add(fld);
    }

    // Ligang: for now copy data from object other, should be interpolation.
    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < size; j++)
        fd(j) = fd_other(j);
    }
  }

// ----------------------------------------------------------------------------

  Increment::Increment(const Increment & other, const bool copy)
    : geom_(new Geometry(*other.geom_)), time_(other.time_),
      vars_(other.vars_) {

    if (vars_.size() != 1) {
      util::abor1_cpp("Increment::Increment(), vars_.size() != 1.",
                        __FILE__, __LINE__);
    }

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<double>(
                         name(var));
      auto fd = make_view<double, 1>(fld);
      fd.assign(0.0);

      atlasFieldSet_->add(fld);
    }

    if (copy) {
      const int size = geom_->atlasFunctionSpace()->size();

      for (int i = 0; i < vars_.size(); i++) {
        auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
        auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));
        for (int j = 0; j < size; j++)
          fd(j) = fd_other(j);
      }
    }
  }

// ----------------------------------------------------------------------------

  Increment::Increment(const Increment & other)
    : geom_(new Geometry(*other.geom_)), time_(other.time_),
      vars_(other.vars_) {

    if (vars_.size() != 1) {
      util::abor1_cpp("Increment::Increment(), vars_.size() != 1.",
                        __FILE__, __LINE__);
    }

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<double>(
                         name(var));
      auto fd = make_view<double, 1>(fld);
      fd.assign(0.0);

      atlasFieldSet_->add(fld);
    }

    // Ligang: copy data from other.
    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < size; j++)
        fd(j) = fd_other(j);
    }
  }
// ----------------------------------------------------------------------------

  Increment::~Increment() { }

// ----------------------------------------------------------------------------

  Increment & Increment::operator =(const Increment &other) {
    // Ligang: Do we assume the geom_ and vars_ are the same, (so they are not
    // empty)? For now yes, as all of the constructors allocate atlasFieldSet.
    // so just assign the values of FieldSet. same for other math operators.
    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < size; j++)
        fd(j) = fd_other(j);
    }

    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator -=(const Increment &other) {
    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < size; j++)
        fd(j) -= fd_other(j);
    }

    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator +=(const Increment &other) {
    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < size; j++)
        fd(j) += fd_other(j);
    }

    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator *=(const double &zz) {
    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
      for (int j = 0; j < size; j++)
        fd(j) *= zz;
    }

    return *this;
  }

// ----------------------------------------------------------------------------

  void Increment::accumul(const double & zz, const State & xx) {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_xx = make_view<double, 1>(xx.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < size; i++)
      fd(i) += zz*fd_xx(i);

    return;
  }

// ----------------------------------------------------------------------------

  void Increment::axpy(const double &zz, const Increment &dx, const bool check)
  {
    ASSERT(!check || time_ == dx.validTime());

    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_dx = make_view<double, 1>(dx.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < size; i++)
      fd(i) += zz*fd_dx(i);

    return;
  }

// ----------------------------------------------------------------------------

  void Increment::diff(const State & x1, const State & x2) {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_x1 = make_view<double, 1>(x1.atlasFieldSet()->field(0));
    auto fd_x2 = make_view<double, 1>(x2.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < size; i++)
      fd(i) = fd_x1(i) - fd_x2(i);

    return;
  }

// ----------------------------------------------------------------------------

  double Increment::dot_product_with(const Increment &other) const {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();
    double dp = 0.0;

    // Ligang: will be updated with missing_value process!
    for (int i = 0; i < size; i++)
      dp += fd(i)*fd_other(i);

    return dp;
  }

// ----------------------------------------------------------------------------

  double Increment::norm() const {
    double norm = 0.0, s = 0.0;

    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    const int size = geom_->atlasFunctionSpace()->size();

    int nValid = 0;
    for (int i = 0; i < size; i++) {
      // Ligang: will be updated with missing_value process!
      if (true) {
        nValid += 1;
        s += fd(i)*fd(i);
      }
    }

    if (nValid == 0)
      norm = 0.0;
    else
      norm = sqrt(s/(1.0*nValid));

    return norm;
  }

// ----------------------------------------------------------------------------

  void Increment::random() {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    const int size = geom_->atlasFunctionSpace()->size();

    util::NormalDistribution<double> x(size, 0, 1.0, 1);
//  util::UniformDistribution<double> x(size, -1.0, 1.0);

    for (int i = 0; i < size; i++)
      fd(i) = x[i];
  }

// ----------------------------------------------------------------------------

  void Increment::schur_product_with(const Increment &rhs ) {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_rhs = make_view<double, 1>(rhs.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();
    for (int i = 0; i < size; i++)
      fd(i) *= fd_rhs(i);

    return;
  }

// ----------------------------------------------------------------------------

  void Increment::zero() {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    fd.assign(0.0);
  }

// ----------------------------------------------------------------------------

  void Increment::ones() {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    fd.assign(1.0);
  }

// ----------------------------------------------------------------------------

  void Increment::zero(const util::DateTime & time) {
    zero();
    time_ = time;
  }

// ----------------------------------------------------------------------------

  void Increment::dirac(const eckit::Configuration & conf) {
    // Ligang: Do NOT know what it is, may ask Jim. Idea of coding from JEDI.
    std::vector<int> ixdir(conf.getIntVector("ixdir"));
    std::vector<int> iydir(conf.getIntVector("iydir"));

    // check
    ASSERT(ixdir.size() > 0 && ixdir.size() == iydir.size());
    const int dir_size = ixdir.size();

    // Ligang: This is where we need the field_data to be 2D;
    // How do we make it 2D? should be related to when create the field.
    const int ny = static_cast<int>(geom_->atlasFunctionSpace()->grid().ny());
    const int nx = static_cast<int>(
      ((atlas::RegularLonLatGrid&)(geom_->atlasFunctionSpace()->grid())).nx() );

    for (int i = 0; i < dir_size; i++)
      ASSERT(ixdir[i] < nx && iydir[i] < ny);

    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    for (int i = 0; i < dir_size; i++) {
      int idx = ixdir[i] + iydir[i]*nx;
      fd(idx) = 1.0;
    }
  }

// ----------------------------------------------------------------------------

  void Increment::print(std::ostream & os) const {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    const int size = geom_->atlasFunctionSpace()->size();

    double mean = 0.0, sum = 0.0,
          min = std::numeric_limits<double>::max(),
          max = std::numeric_limits<double>::min();

    // Does fd here (Atlas) provide max, min, mean functions as c++ std does,
    // Since normally Increment does not have missing_values?
    for (int i = 0; i < size; i++) {
      if (fd(i) < min) min = fd(i);
      if (fd(i) > max) max = fd(i);
      sum += fd(i);
    }
    mean = sum / (1.0*size);

    os << "Increment: print diagnostic info about the increment here"
       << "(min/max/mean for var sea_surface_temperature: )"
       << std::endl;
    os << "min = " << min << ", max = " << max << ", mean = " << mean
       << std::endl;
  }

// ----------------------------------------------------------------------------

  void Increment::read(const eckit::Configuration & conf) {
    int time = 0, lon = 0, lat = 0, iread = 0;
    std::string sdate, filename, record;

    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    const int size = geom_->atlasFunctionSpace()->size();

    // get filename
    if (!conf.get("filename", filename))
      util::abor1_cpp("Increment::read(), Get filename failed.",
        __FILE__, __LINE__);

    // open netCDF file
    netCDF::NcFile file(filename.c_str(), netCDF::NcFile::read);
    if (file.isNull())
      util::abor1_cpp("Increment::read(), Create netCDF file failed.",
        __FILE__, __LINE__);

    // get file dimensions
    time = static_cast<int>(file.getDim("time").getSize());
    lon  = static_cast<int>(file.getDim("lon").getSize());
    lat  = static_cast<int>(file.getDim("lat").getSize());
    if (time != 1 ||
        lat != static_cast<int>(geom_->atlasFunctionSpace()->grid().ny()) ||
        lon != static_cast<int>((((atlas::RegularLonLatGrid&)  // LC: no &?
               (geom_->atlasFunctionSpace()->grid()))).nx()) ) {
      util::abor1_cpp("Increment::read(), lat!=ny or lon!=nx",
        __FILE__, __LINE__);
    }

    // get sst increment data
    netCDF::NcVar incVar;
    incVar = file.getVar("sst");  // Ligang: will change the var name!
    if (incVar.isNull())
      util::abor1_cpp("Get inc var failed.", __FILE__, __LINE__);

    double incData[lat][lon];
    incVar.getVar(incData);

    int idx = 0;
    for (int j = 0; j < lat; j++)
      for (int i = 0; i < lon; i++)
        fd(idx++) = incData[j][i];
  }

// ----------------------------------------------------------------------------

  void Increment::write(const eckit::Configuration & conf) const {
    int lat, lon, time = 1;
    std::string filename;

    // Ligang: Just for test, should NOT be here, const! put it outside.
    // compile error, discards qualifiers [-fpermissive], "const"!
    // random();

    // get filename
    if (!conf.get("filename", filename))
      util::abor1_cpp("Increment::write(), Get filename failed.",
                      __FILE__, __LINE__);
    else
      oops::Log::info() << "Increment::write(), filename=" << filename
                        << std::endl;

    // create netCDF file
    netCDF::NcFile file(filename.c_str(), netCDF::NcFile::replace);
    if (file.isNull())
      util::abor1_cpp("Increment::write(), Create netCDF file failed.",
                      __FILE__, __LINE__);

    // define dims
    lat = geom_->atlasFunctionSpace()->grid().ny();
    lon =((atlas::RegularLonLatGrid)(geom_->atlasFunctionSpace()->grid())).nx();

    // unlimited dim if without size parameter, then it'll be 0,
    // what about the size?
    netCDF::NcDim timeDim = file.addDim("time", 1);
    netCDF::NcDim latDim  = file.addDim("lat" , lat);
    netCDF::NcDim lonDim  = file.addDim("lon" , lon);
    if (timeDim.isNull() || latDim.isNull() || lonDim.isNull())
      util::abor1_cpp("Increment::write(), Define dims failed.",
                      __FILE__, __LINE__);

    std::vector<netCDF::NcDim> dims;
    dims.push_back(timeDim);
    dims.push_back(latDim);
    dims.push_back(lonDim);

    // Lignag: define coordinate vars "lat" and "lon"

    // Ligang: define units atts for coordinate vars

    // define data vars, will change the name later.

    netCDF::NcVar incVar = file.addVar(std::string("increment")
                         , netCDF::ncDouble, dims);

    // Ligang: define units atts for data vars
    incVar.putAtt("units", "K");
    incVar.putAtt("_FillValue", netCDF::NcDouble(), -32768.);
    incVar.putAtt("missing_value", netCDF::NcDouble(), -32768.);


    // write data to the file
    // Ligang: compile failed, rank should be 1, not 2, related to init?
//  auto field_data = make_view<double, 2>(atlasFieldSet_->field(0));
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));

    double incData[time][lat][lon];
    int idx = 0;
    for (int i = 0; i < lat; i++)
      for (int j = 0; j < lon; j++)
        incData[0][i][j] = fd(idx++);

    incVar.putVar(incData);

    oops::Log::info() << "Increment::write(), Successfully write data to file!"
                      << std::endl;
  }

// ----------------------------------------------------------------------------
}  // namespace oisst
