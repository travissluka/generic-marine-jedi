/*
 * (C) Copyright 2021-2021 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "umdsst/VariableChange/Model2GeoVaLs.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/util/abor1_cpp.h"

namespace umdsst {

// -----------------------------------------------------------------------------

static oops::VariableChangeMaker<Traits, Model2GeoVaLs>
  makerVariableChangeModel2GeoVaLs_("Model2GeoVaLs");

static oops::VariableChangeMaker< Traits, Model2GeoVaLs>
  makerVariableChangeDefault_("default");

// -----------------------------------------------------------------------------

Model2GeoVaLs::Model2GeoVaLs(const Geometry & geom,
                             const eckit::Configuration & conf)
  : geom_(new Geometry(geom)) {
}

// -----------------------------------------------------------------------------

void Model2GeoVaLs::changeVar(const State & xin, State & xout) const {
  const int size = geom_->atlasFunctionSpace()->size();

  for ( int i = 0; i < xout.variables().size(); i++ ) {
    std::string name = xout.variables()[i];

    if ( xin.variables().has(name) ) {
      // variable is simply being copied from xin to xout
      xout.atlasFieldSet()->field(name) = xin.atlasFieldSet()->field(name);

    } else if (name == "sea_area_fraction") {
      // convert integer land mask to a floating point field
      auto fd = atlas::array::make_view<double, 2>(
          xout.atlasFieldSet()->field(name));
      auto fd_src = atlas::array::make_view<int, 2>(
          geom_->atlasFieldSet()->field("gmask"));
      for (int j=0; j < size; j++)
         fd(j, 0) = static_cast<double>(fd_src(j, 0));
    }
  }
}

// -----------------------------------------------------------------------------

void Model2GeoVaLs::changeVarInverse(const State &, State &) const {
  util::abor1_cpp("Model2GeoVaLs::changeVarInverse not implemented");
}

// -----------------------------------------------------------------------------

}  // namespace umdsst
