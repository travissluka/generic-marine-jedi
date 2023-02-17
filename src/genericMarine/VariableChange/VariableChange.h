/*
 * (C) Copyright 2021-2022 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/base/VariableChangeParametersBase.h"

namespace genericMarine {

  // Forward declarations
  class Geometry;
  class State;

// -----------------------------------------------------------------------------

class VariableChangeParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters,
                           oops::VariableChangeParametersBase)
 public:
  // Wrapper to VariableChange parameters
//   VariableChangeParametersWrapper variableChangeParametersWrapper{this};
};

// -----------------------------------------------------------------------------

class VariableChange : public util::Printable {
 public:
  static const std::string classname() {return "genericMarine::VariableChange";}
  typedef VariableChangeParameters Parameters_;

  explicit VariableChange(const Parameters_ &, const Geometry &);
  ~VariableChange() = default;

  void changeVar(State &, const oops::Variables &) const;
  void changeVarInverse(const State &, const oops::Variables &) const { ASSERT(1 == 2); }

 private:
  void print(std::ostream &) const override { ASSERT(1 == 2); }
  const Geometry & geom_;
};

}  // namespace genericMarine
