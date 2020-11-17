/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "umdsst/ModelAux/ModelAuxControl.h"
#include "umdsst/ModelAux/ModelAuxIncrement.h"

#include "oops/util/abor1_cpp.h"

namespace umdsst {

// ----------------------------------------------------------------------------

  ModelAuxControl::ModelAuxControl(const Geometry &,
                                   const eckit::Configuration &) {
  }

// ----------------------------------------------------------------------------

  ModelAuxControl::ModelAuxControl(const Geometry &,
                                   const ModelAuxControl &) {
  }

// ----------------------------------------------------------------------------

  ModelAuxControl::ModelAuxControl(const ModelAuxControl &, const bool) {
  }

// ----------------------------------------------------------------------------

  ModelAuxControl::~ModelAuxControl() {
  }

// ----------------------------------------------------------------------------

  ModelAuxControl & ModelAuxControl::operator+=(const ModelAuxIncrement & inc) {
    return *this;
  }
// ----------------------------------------------------------------------------

  void ModelAuxControl::print(std::ostream & os) const {
    os << "umdsst::ModelAuxControl not implemented"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace umdsst
