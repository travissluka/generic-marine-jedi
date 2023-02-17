/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Geometry/Geometry.h"
#include "genericMarine/ModelAux/ModelAuxCovariance.h"

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"

namespace genericMarine {

// ----------------------------------------------------------------------------

  ModelAuxCovariance::ModelAuxCovariance(const eckit::Configuration & conf,
                                         const Geometry & geom)
    : conf_(conf) {
  }

// ----------------------------------------------------------------------------

  ModelAuxCovariance::~ModelAuxCovariance() {
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::linearize(const ModelAuxControl &,
                                     const Geometry &) {
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::multiply(const ModelAuxIncrement &,
                                    ModelAuxIncrement &) {
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::inverseMultiply(const ModelAuxIncrement &,
                                           ModelAuxIncrement &) const {
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::randomize(ModelAuxIncrement &) const {
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::print(std::ostream & os) const {
    os << "genericMarine::ModelAuxCovariance not implemented"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace genericMarine
