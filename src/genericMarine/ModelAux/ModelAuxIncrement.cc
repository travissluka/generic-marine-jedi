/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/ModelAux/ModelAuxIncrement.h"
#include "genericMarine/ModelAux/ModelAuxControl.h"

#include "oops/util/abor1_cpp.h"

namespace genericMarine {

// ----------------------------------------------------------------------------

  ModelAuxIncrement::ModelAuxIncrement(const ModelAuxIncrement &,
                                       const eckit::Configuration &) {
  }

// ----------------------------------------------------------------------------

  ModelAuxIncrement::ModelAuxIncrement(const ModelAuxIncrement &, const bool) {
  }

// ----------------------------------------------------------------------------

  ModelAuxIncrement::ModelAuxIncrement(const Geometry &,
                                       const eckit::Configuration &) {
  }

// ----------------------------------------------------------------------------

  ModelAuxIncrement::~ModelAuxIncrement() {
  }

// ----------------------------------------------------------------------------
  void ModelAuxIncrement::diff(const ModelAuxControl &,
                               const ModelAuxControl &) {
  }
// ----------------------------------------------------------------------------

  void ModelAuxIncrement::zero() {
  }

// ----------------------------------------------------------------------------

  ModelAuxIncrement & ModelAuxIncrement::operator*=(const double) {
    return *this;
  }

// ----------------------------------------------------------------------------

  ModelAuxIncrement & ModelAuxIncrement::operator+=(const ModelAuxIncrement &) {
    return *this;
  }

// ----------------------------------------------------------------------------

  ModelAuxIncrement & ModelAuxIncrement::operator-=(const ModelAuxIncrement &) {
    return *this;
  }

// ----------------------------------------------------------------------------

  double ModelAuxIncrement::norm() const {
    return 0.0;
  }

// ----------------------------------------------------------------------------

  void ModelAuxIncrement::axpy(const double, const ModelAuxIncrement &) {
  }

// ----------------------------------------------------------------------------

  double ModelAuxIncrement::dot_product_with(const ModelAuxIncrement &) const {
    return 0.0;
  }

// ----------------------------------------------------------------------------

  void ModelAuxIncrement::print(std::ostream & os) const {
    os << "genericMarine::ModelAuxIncrement not implemented"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace genericMarine
