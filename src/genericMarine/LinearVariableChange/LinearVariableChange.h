/*
 * (C) Copyright 2022-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/util/Printable.h"

// forward declartions
namespace oops {
  class Variables;
}
namespace genericMarine {
  // forward declarations
  class State;
  class Geometry;
  class Increment;
}

// -----------------------------------------------------------------------------

namespace genericMarine {

// -----------------------------------------------------------------------------

class LinearVariableChangeParameters : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters,
                           oops::LinearVariableChangeParametersBase)
};

// -----------------------------------------------------------------------------

class LinearVariableChange : public util::Printable {
  // NOTE: this is an empty class that doesn't do anything
 public:
  static const std::string classname() {return "genericMarine::LinearVariableChange";}
  typedef LinearVariableChangeParameters Parameters_;

  explicit LinearVariableChange(const Geometry &, const Parameters_ &) {}
  ~LinearVariableChange() {}

  void changeVarTraj(const State &, const oops::Variables &) {}

  void changeVarTL(Increment &, const oops::Variables &) const {}
  void changeVarInverseTL(Increment &, const oops::Variables &) const {}
  void changeVarAD(Increment &, const oops::Variables &) const {}
  void changeVarInverseAD(Increment &, const oops::Variables &) const {}

 private:
  void print(std::ostream &) const override { ASSERT(1 == 2); }
};

}  // namespace genericMarine
