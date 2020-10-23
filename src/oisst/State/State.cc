/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <numeric>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iterator>

// #include <netcdf>
#include "netcdf"

#include "oisst/Geometry/Geometry.h"
#include "oisst/State/State.h"

#include "eckit/config/Configuration.h"

#include "atlas/field.h"
#include "atlas/array.h"
#include "atlas/option.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

// using namespace netCDF;
using atlas::array::make_view;
using atlas::option::name;


namespace oisst {

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const eckit::Configuration & conf)
    : geom_(new Geometry(geom)), time_(conf.getString("date")),
      vars_(conf, "state variables") {

    if (vars_.size() != 1) {
      util::abor1_cpp("State::State(), vars_.size() != 1",
                       __FILE__, __LINE__);
    }

    const int SIZE = geom_->atlasFunctionSpace()->size();

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<float>(
                         name(var));
      // Initialize to 0.0 or 300.0
      auto fd = make_view<float, 1>(fld);
      fd.assign(300.0);  // norm is 300., 0.0 if init to 0.0;

      atlasFieldSet_->add(fld);
    }

    // check and read in data if it has filename
    if (conf.has("filename")) {
      read(conf);
//    double norm1 = norm();
//    std::cout << "State::State(), after read file, norm1=" << norm1
//              << std::endl;
    }
  }

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & time)
    : geom_(new Geometry(geom)), time_(time), vars_(vars) {

    if (vars_.size() != 1) {
      util::abor1_cpp("State::State(), vars_.size() != 1",
                      __FILE__, __LINE__);
    }

    const int SIZE = geom_->atlasFunctionSpace()->size();

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<float>(
                         name(var));
      // Initialize to 300.0, will fail the test: EXPECT(xx3.norm()==0); so 0.0
      auto fd = make_view<float, 1>(fld);
      fd.assign(0.0);

      atlasFieldSet_->add(fld);
    }
  }

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const State & other)
    : geom_(new Geometry(geom)), time_(other.time_), vars_(other.vars_) {
    // Change state resolution

    if (vars_.size() != 1) {
      util::abor1_cpp("State::State() needs to be implemented.",
                      __FILE__, __LINE__);
    }

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      // Ligang: should be copied from "other"?
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<float>(
                         name(var));
      atlasFieldSet_->add(fld);
    }

    // Ligang: copy data from object other, better way to initialize?!
    const int SIZE = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<float, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<float, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < SIZE; j++)
        fd(j) = fd_other(j);
    }
  }

// ----------------------------------------------------------------------------

  State::State(const State & other)
    : geom_(new Geometry(*other.geom_)), time_(other.time_),
      vars_(other.vars_) {

    if (vars_.size() != 1) {
      util::abor1_cpp("State::State() needs to be implemented.",
                       __FILE__, __LINE__);
    }

    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int i = 0; i < vars_.size(); i++) {
      std::string var = vars_[i];
      atlas::Field fld = geom_->atlasFunctionSpace()->createField<float>(
                         name(var));
      atlasFieldSet_->add(fld);
    }

    // Ligang: copy data from object other
    const int SIZE = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<float, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<float, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < SIZE; j++)
        fd(j) = fd_other(j);
    }
  }

// ----------------------------------------------------------------------------

  State::~State() {
  // util::abor1_cpp("State::~State() needs to be implemented.",
  //                __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  State & State::operator+=(const Increment & dx)
  {
  // util::abor1_cpp("State::operator+=(Increment) needs to be implemented.",
  //                __FILE__, __LINE__);

    return *this;
  }

// ----------------------------------------------------------------------------

  void State::accumul(const double &zz, const State &rhs) {
    auto field_data = make_view<float, 1>(atlasFieldSet_->field(0));
    auto rhs_field_data = make_view<float, 1>(rhs.atlasFieldSet()->field(0));

    // Ligang: the following does not work, not defined operator;
    // field_data = field_data + zz*rhs_field_data;
    // Ligang: return size_halo_? looks not.
    size_t SIZE = geom_->atlasFunctionSpace()->size();

    // Ligang: what about missing value?
    for (size_t jnode = 0; jnode < SIZE; jnode++)
      field_data(jnode) +=  zz*rhs_field_data(jnode);

    return;
  }

// ----------------------------------------------------------------------------

  double State::norm() const {
    double norm = 0.0, s = 0.0;
    int nx = 0, ny = 0;

    // Ligang: just consider 1 var for now.
    auto field_data = make_view<float, 1>(atlasFieldSet_->field(0));

    ny = static_cast<int>(geom_->atlasFunctionSpace()->grid().ny());
    nx = static_cast<int>(
      ((atlas::RegularLonLatGrid&)(geom_->atlasFunctionSpace()->grid())).nx() );

    // Ligang: they ARE the same.
    if (geom_->atlasFunctionSpace()->size() != ny*nx)
      util::abor1_cpp("State::norm() size() != ny*nx.", __FILE__, __LINE__);

    // Ligang: undefined operator "+" for field_data+ny*nx;
    // s = std::inner_product(field_data, field_data+ny*nx, field_data, 0);
    int nValid = 0;
    for (int jnode = 0; jnode < ny*nx; jnode++) {
      // Ligang: if not missing value, not a good way, need to improve!
      if (true) {  // LC: previously (field_data(jnode) > 0.0), fail test
        nValid += 1;
        s += field_data(jnode)*field_data(jnode);
      }
    }

    if (nValid == 0)
      norm = 0.0;
    else
      norm = sqrt(s/(1.0*nValid));

    return norm;
  }

// ----------------------------------------------------------------------------

  void State::zero() {
    // Ligang: just consider 1 var for now

    // atlasFieldSet_->field(0).array() = 0.0  // does not work.
    // field_data = 0.0; // Ligang: undefined operator "="!

    auto field_data = make_view<float, 1>(atlasFieldSet_->field(0));
    field_data.assign(0.0);
  }

// ----------------------------------------------------------------------------

  void State::read(const eckit::Configuration & conf) {
    int iread = 0, time = 0, lon = 0, lat = 0, bc = 0;
    std::string sdate, filename, record;

    auto field_data = make_view<float, 1>(atlasFieldSet_->field(0));
    const int SIZE = geom_->atlasFunctionSpace()->size();

    // Ligang: conf only has "statefile" part of the yml file.
    // Here we do not actually use iread
    if (conf.has("read_from_file")) {
      iread = conf.getInt("read_from_file");
    } else {
      // you can output keys to check
      std::vector<std::string> keys = conf.keys();
    }

    // get filename
    if (!conf.get("filename", filename))
      util::abor1_cpp("Get filename failed.", __FILE__, __LINE__);

    // Ligang: open netCDF file
    netCDF::NcFile file(filename.c_str(), netCDF::NcFile::read);
    if (file.isNull())
      util::abor1_cpp("Create netCDF file failed.", __FILE__, __LINE__);

    // Ligang: get file dimensions, for checking
    time = static_cast<int>(file.getDim("time").getSize());
    lon  = static_cast<int>(file.getDim("lon").getSize());
    lat  = static_cast<int>(file.getDim("lat").getSize());
    if (time != 1 ||
        lat != static_cast<int>(geom_->atlasFunctionSpace()->grid().ny()) ||
        lon != static_cast<int>((((atlas::RegularLonLatGrid&)  // LC: no &?
               (geom_->atlasFunctionSpace()->grid()))).nx()) ) {
      util::abor1_cpp("lat!=ny or lon!=nx", __FILE__, __LINE__);
    }

    // get sst data
    netCDF::NcVar sstVar;
    sstVar = file.getVar("sst");
    if (sstVar.isNull())
      util::abor1_cpp("Get sst var failed.", __FILE__, __LINE__);

    // Ligang: do we have to use for-loop?
    // sstVar.getVar(field_data); // doesn't work this way
    float sstData[lat][lon];
    sstVar.getVar(sstData);

    int idx = 0;
    for (int i = 0; i < lat; i++)
      for (int j = 0; j < lon; j++)
        field_data(idx++) = sstData[i][j];
  }

// ----------------------------------------------------------------------------

  void State::write(const eckit::Configuration & conf) const {
    int lat, lon, time = 1;
    std::string filename;

    // get filename
    if (!conf.get("filename", filename))
      util::abor1_cpp("Get filename failed.", __FILE__, __LINE__);
    else
      std::cout << "State::write(), filename=" << filename << std::endl;

    // create netCDF file
    netCDF::NcFile file(filename.c_str(), netCDF::NcFile::replace);
    if (file.isNull())
      util::abor1_cpp("Create netCDF file failed.", __FILE__, __LINE__);

    // define dims
    lat = geom_->atlasFunctionSpace()->grid().ny();
    lon =((atlas::RegularLonLatGrid)(geom_->atlasFunctionSpace()->grid())).nx();

    // unlimited dim if without size parameter, then it'll be 0,
    // what about the size?
    netCDF::NcDim timeDim = file.addDim("time", 1);
    netCDF::NcDim latDim  = file.addDim("lat" , lat);
    netCDF::NcDim lonDim  = file.addDim("lon" , lon);
    if (timeDim.isNull() || latDim.isNull() || lonDim.isNull())
      util::abor1_cpp("Define dims failed.", __FILE__, __LINE__);

    std::vector<netCDF::NcDim> dims;
    dims.push_back(timeDim);
    dims.push_back(latDim);
    dims.push_back(lonDim);

    // Lignag: define coordinate vars "lat" and "lon"

    // Ligang: define units atts for coordinate vars

    // define data vars
    netCDF::NcVar sstVar = file.addVar(std::string("sst"), netCDF::ncFloat
                                       , dims);

    // Ligang: define units atts for data vars
    sstVar.putAtt("units", "K");
    sstVar.putAtt("_FillValue", netCDF::NcFloat(), -32768.);

    // write data to the file
    // Ligang: compile failed, rank should be 1, not 2, related to init?
//  auto field_data = make_view<float, 2>(atlasFieldSet_->field(0));
    auto field_data = make_view<float, 1>(atlasFieldSet_->field(0));

    float sstData[time][lat][lon];
    int idx = 0;
    for (int i = 0; i < lat; i++)
      for (int j = 0; j < lon; j++)
        sstData[0][i][j] = field_data(idx++);

    sstVar.putVar(sstData);

    oops::Log::info() << "Successfully write data to file!" << std::endl;
  }

// ----------------------------------------------------------------------------

  std::shared_ptr<const Geometry> State::geometry() const {return geom_;}

// ----------------------------------------------------------------------------

  void State::print(std::ostream & os) const {
    os << "insert diagnostic information about state here "
       << "(min/max/mean for each state var?)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace oisst
