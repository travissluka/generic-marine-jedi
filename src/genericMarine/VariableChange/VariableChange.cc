/*
 * (C) Copyright 2021-2022 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Geometry/Geometry.h"
#include "genericMarine/State/State.h"
#include "genericMarine/VariableChange/VariableChange.h"

#include "oops/util/abor1_cpp.h"

namespace genericMarine {

// -----------------------------------------------------------------------------

VariableChange::VariableChange(const Parameters_ & params,
                               const Geometry & geom)
    : geom_(geom) {
}

// -----------------------------------------------------------------------------

void VariableChange::changeVar(State & xx, const oops::Variables & vars) const {
  State xx2(xx.geometry(), vars, xx.time());

  for ( int v = 0; v < vars.size(); v++ ) {
    std::string name = vars[v];
    const int size = geom_.functionSpace().size();

    auto dst = atlas::array::make_view<double, 2>(xx2.fieldSet().field(name));

    if ( xx.variables().has(name) ) {
      // no variable change required, already exists
      auto src = atlas::array::make_view<double, 2>(xx.fieldSet().field(name));
      for (int j=0; j < size; j++) dst(j, 0) = src(j, 0);

    } else if (name == "sea_area_fraction") {
      // convert integer land mask to a floating point field
      auto src = atlas::array::make_view<int, 2>(
        geom_.extraFields().field("gmask"));
      for (int j=0; j < size; j++) {
        dst(j, 0) = static_cast<double>(src(j, 0));
      }
    } else {
      ASSERT(1 == 2);
    }

    xx2.fieldSet().field(name).metadata().set("interp_type", "default");
  }

  xx = xx2;
}

// -----------------------------------------------------------------------------

}  // namespace genericMarine
