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

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace umdsst {

// ----------------------------------------------------------------------------

  Covariance::Covariance(const Geometry & geom,
                         const oops::Variables & vars,
                         const eckit::Configuration & conf,
                         const State & x1, const State & x2) {
    oops::Log::trace() << "Covariance created!" << std::endl;  // ?
  }

// ----------------------------------------------------------------------------

  Covariance::~Covariance() {
    oops::Log::trace() << "Covariance destructed!" << std::endl;
  }

// ----------------------------------------------------------------------------

  void Covariance::inverseMultiply(const Increment & dxin,
                                   Increment & dxout) const {
    dxout = dxin;  // Identity Cov
  }

// ----------------------------------------------------------------------------

  void Covariance::multiply(const Increment & dxin, Increment & dxout) const {
    dxout = dxin;  // Identity
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
