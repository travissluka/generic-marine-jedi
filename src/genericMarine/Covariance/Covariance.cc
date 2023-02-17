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


namespace genericMarine {

// ----------------------------------------------------------------------------

  Covariance::Covariance(const Geometry & geom,
                         const oops::Variables & vars,
                         const eckit::Configuration & conf,
                         const State & x1, const State & x2) {
    oops::Log::trace() << "genericMarine::Covariance::Covariance starting"<< std::endl;
  }

// ----------------------------------------------------------------------------

  Covariance::~Covariance() {
    oops::Log::trace() << "Covariance destructed!" << std::endl;
  }

// ----------------------------------------------------------------------------

  void Covariance::inverseMultiply(const Increment & dxin, Increment & dxout) const {
    dxout = dxin;
  }

// ----------------------------------------------------------------------------

  void Covariance::multiply(const Increment & dxin, Increment & dxout) const {
    dxout = dxin;
  }

// ----------------------------------------------------------------------------

  void Covariance::randomize(Increment & dx) const {
    dx.random();
  }

// ----------------------------------------------------------------------------

  void Covariance::print(std::ostream & os) const {
    os << "Covariance::print not implemented" << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace genericMarine
