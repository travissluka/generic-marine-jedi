/*
 * (C) Copyright 2021-2021 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/VariableChangeBase.h"
#include "oops/generic/IdVariableChange.h"
#include "umdsst/Traits.h"

namespace umdsst {

static oops::GenericVariableChangeMaker<Traits,
  oops::IdVariableChange<Traits> > makerVariableChange_("default");

}  // namespace umdsst
