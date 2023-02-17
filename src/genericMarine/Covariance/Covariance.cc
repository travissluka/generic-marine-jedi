/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <limits>
#include <ostream>
#include <string>

#include "genericMarine/Covariance/Covariance.h"
#include "genericMarine/Geometry/Geometry.h"
#include "genericMarine/Increment/Increment.h"
#include "genericMarine/State/State.h"

#include "eckit/config/Configuration.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "saber/bump/type_bump.h"

namespace genericMarine {

// ----------------------------------------------------------------------------

  Covariance::Covariance(const Geometry & geom,
                         const oops::Variables & vars,
                         const eckit::Configuration & conf,
                         const State & x1, const State & x2) {
    oops::Log::trace() << "genericMarine::Covariance::Covariance starting"<< std::endl;

    // setup BUMP
    // --------------------------------------------
    eckit::LocalConfiguration bumpConf;
    conf.get("bump", bumpConf);
    const double msvalr = util::missingValue(msvalr);
    bumpConf.set("msvalr", msvalr);

    eckit::LocalConfiguration gridConf;
    std::string prefix;
    bumpConf.get("prefix", prefix);
    gridConf.set("prefix", prefix + "_00");
    gridConf.set("variables", vars.variables());
    gridConf.set("nv", vars.size());
    gridConf.set("nl0", 1);

    oops::Log::info() << "Configuration: " << bumpConf << std::endl;
    oops::Log::info() << "Grid " << 0 << ": " << gridConf << std::endl;

    atlas::FieldSet universe_rad = atlas::FieldSet();

    saber::bump_create_f90(keyBump_, &geom.getComm(),
                    geom.functionSpace().get(),
                    geom.extraFields().get(),
                    bumpConf, gridConf, universe_rad.get());

    // pass user generated fields to BUMP
    // NOTE: assuming only 1 state variable
    // --------------------------------------------
    assert(vars.size() == 1);
    std::string param_name;
    atlas::FieldSet param_fieldSet;
    atlas::Field param_field = geom.functionSpace().createField<double>(
      atlas::option::levels(1) | atlas::option::name(vars[0]));
    param_fieldSet.add(param_field);
    auto param_view = atlas::array::make_view<double, 2>(param_field);

    // horizontal correlation lengths
    if (conf.has("correlation lengths")) {
      eckit::LocalConfiguration corrConf;
      conf.get("correlation lengths", corrConf);

      // rh is calculated as follows :
      // 1) rh = "base value" + rossby_radius * "rossby mult"
      // 2) minimum value of "min grid mult" * grid_size is imposed
      // 3) min/max are imposed based on "min value" and "max value"
      // 4) converted from a gaussian sigma to Gaspari-Cohn cutoff distance
      double baseValue = corrConf.getDouble("base value", 0.0);
      double rossbyMult = corrConf.getDouble("rossby mult", 0.0);
      double minGridMult = corrConf.getDouble("min grid mult", 0.0);
      double minValue = corrConf.getDouble("min value", 0.0);
      double maxValue = corrConf.getDouble("max value",
                                       std::numeric_limits<double>::max());

      auto rossbyRadius = atlas::array::make_view<double, 2>(
        geom.extraFields().field("rossby_radius"));
      auto area = atlas::array::make_view<double, 2>(
         geom.extraFields().field("area"));

      param_view.assign(baseValue);
      for ( int i = 0; i < param_field.size(); i++ ) {
        param_view(i, 0) += rossbyMult * rossbyRadius(i, 0);
        param_view(i, 0) = std::max(param_view(i, 0),
                                    minGridMult*sqrt(area(i, 0)));
        param_view(i, 0) = std::max(param_view(i, 0), minValue);
        param_view(i, 0) = std::min(param_view(i, 0), maxValue);
      }

      // note: BUMP expects the length as a Gaspari-Cohn cutoff length,
      //   but we probably think of it as a Gaussian 1 sigma, so convert.
      for ( int i = 0; i < param_field.size(); i++ ) {
        param_view(i, 0) *= 3.57;  // gaussian to GC factor
      }

      param_name = "rh";
      saber::bump_set_parameter_f90(keyBump_, param_name.size(),
                                    param_name.c_str(), 1, param_fieldSet.get());

      // vertical lengths (leave at 1.0, because we have no vertical)
      param_name = "rv";
      param_view.assign(1.0);
      saber::bump_set_parameter_f90(keyBump_, param_name.size(),
                                    param_name.c_str(), 1, param_fieldSet.get());
    }

    // Calculate static B and cleanup
    // --------------------------------------------
    saber::bump_run_drivers_f90(keyBump_);
    saber::bump_partial_dealloc_f90(keyBump_);

    oops::Log::trace() << "genericMarine::Covariance::Covariance done" << std::endl;
  }

// ----------------------------------------------------------------------------

  Covariance::~Covariance() {
    oops::Log::trace() << "Covariance destructed!" << std::endl;
  }

// ----------------------------------------------------------------------------

  void Covariance::inverseMultiply(const Increment & dxin,
                                   Increment & dxout) const {
    oops::IdentityMatrix<Increment> Id;
    dxout.zero();
    GMRESR(dxout, dxin, *this, Id, 10, 1.0e-3);
  }

// ----------------------------------------------------------------------------

  void Covariance::multiply(const Increment & dxin, Increment & dxout) const {
    dxout = dxin;
    saber::bump_apply_nicas_f90(keyBump_, dxout.fieldSet().get());
  }

// ----------------------------------------------------------------------------

  void Covariance::randomize(Increment & dx) const {
    dx.random();
  }

// ----------------------------------------------------------------------------

  void Covariance::print(std::ostream & os) const {
    os << "(TODO, print diagnostic info about the Covariance here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace genericMarine
