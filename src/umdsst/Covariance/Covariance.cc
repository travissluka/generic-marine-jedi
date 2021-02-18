/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "umdsst/Covariance/Covariance.h"
#include "umdsst/Geometry/Geometry.h"
#include "umdsst/Increment/Increment.h"
#include "umdsst/State/State.h"

#include "eckit/config/Configuration.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "saber/bump/type_bump.h"

namespace umdsst {

// ----------------------------------------------------------------------------

  Covariance::Covariance(const Geometry & geom,
                         const oops::Variables & vars,
                         const eckit::Configuration & conf,
                         const State & x1, const State & x2) {
    oops::Log::trace() << "umdsst::Covariance::Covariance starting"<< std::endl;

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
    gridConf.set("nl", 1);

    oops::Log::info() << "Configuration: " << bumpConf << std::endl;
    oops::Log::info() << "Grid " << 0 << ": " << gridConf << std::endl;

    saber::bump_create_f90(keyBump_, &geom.getComm(),
                    geom.atlasFunctionSpace()->get(),
                    geom.atlasFieldSet()->get(),
                    bumpConf, gridConf);

    // pass user generated fields to BUMP
    // NOTE: assuming only 1 state variable
    // --------------------------------------------
    assert(vars.size() == 1);
    std::string param_name;
    atlas::FieldSet param_fieldSet;
    atlas::Field param_field = geom.atlasFunctionSpace()->createField<double>(
      atlas::option::levels(1) | atlas::option::name(vars[0]));
    param_fieldSet.add(param_field);
    auto param_view = atlas::array::make_view<double, 2>(param_field);

    // horizontal correlation lengths
    if (conf.has("correlation lengths")) {
      eckit::LocalConfiguration corrConf;
      conf.get("correlation lengths", corrConf);
      if ( corrConf.has("fixed") ) {
        // single global value for correlation length given
        double val;
        corrConf.get("fixed", val);
        param_view.assign(val);
      } else {
        util::abor1_cpp("Covariance::Covariance() no correlation length "
                        "method specified", __FILE__, __LINE__);
      }
      param_name = "cor_rh";
      saber::bump_set_parameter_f90(keyBump_, param_name.size(),
                                    param_name.c_str(), param_fieldSet.get());

      // vertical lengths (leave at 1.0, because we have no vertical)
      param_name = "cor_rv";
      param_view.assign(1.0);
      saber::bump_set_parameter_f90(keyBump_, param_name.size(),
                                    param_name.c_str(), param_fieldSet.get());
    }

    // Calculate static B and cleanup
    // --------------------------------------------
    saber::bump_run_drivers_f90(keyBump_);
    saber::bump_partial_dealloc_f90(keyBump_);

    oops::Log::trace() << "umdsst::Covariance::Covariance done" << std::endl;
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
    saber::bump_apply_nicas_f90(keyBump_, (*dxout.atlasFieldSet()).get());
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

}  // namespace umdsst
