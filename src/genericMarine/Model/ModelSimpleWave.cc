/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Traits.h"
#include "genericMarine/Model/ModelSimpleWave.h"
#include "genericMarine/Geometry/Geometry.h"

#include "oops/interface/ModelBase.h"
#include "oops/util/FieldSetHelpers.h"

namespace genericMarine {

  static oops::interface::ModelMaker<Traits, ModelSimpleWave> modelSimpleWave_("SimpleWave");

// -----------------------------------------------------------------------------

  ModelSimpleWave::ModelSimpleWave(const Geometry & geom, const ModelSimpleWaveParameters & params)
  : geom_(geom), tstep_(params.tstep) {

    // calculate required derived fields
    dx_ = geom.functionSpace().createField<double>(
      atlas::option::levels(1) | atlas::option::name("dx"));
    auto fd = atlas::array::make_view<double, 2>(dx_);
    fd.assign(110e3);  //TODO replace this with actual value
  }

// -----------------------------------------------------------------------------

  void ModelSimpleWave::print(std::ostream & os) const {
    os << "ModelSimpleWave::print not implemented";
  }

// -----------------------------------------------------------------------------

  void ModelSimpleWave::initialize(State & xx) const {
    fs_m1_.clear();
  }

// -----------------------------------------------------------------------------

void ModelSimpleWave::step(State & xx, const ModelAuxControl &) const {
  // grid related stuff
  atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  atlas::StructuredGrid grid = fspace.grid();
  const int size = grid.size();
  atlas::idx_t ny = grid.ny();

  // copy the existing data
  atlas::FieldSet fs_p1 = xx.fieldSet();
  atlas::FieldSet fs = util::copyFieldSet(fs_p1);

  // get views to the fields we'll need
  auto c_view = atlas::array::make_view<double,2>(geom_.extraFields().field("phase_speed"));
  auto halo_view = atlas::array::make_view<int, 1>(fspace.ghost());

  for (int v = 0; v < fs.size(); v++) {
    auto v_view = atlas::array::make_view<double,2>(fs[v]);
    auto v_p1_view = atlas::array::make_view<double,2>(fs_p1[v]);

    geom_.functionSpace().haloExchange(fs[v]);
    geom_.functionSpace().haloExchange(fs_p1[v]);
    // WRONG
    // for(atlas::idx_t idx = 0; idx < v_view.size(); idx++) {
    //   if (halo_view(idx))continue;

    //   atlas::idx_t x, y;
    //   grid.index2ij(idx, x, y);
    //   atlas::idx_t idx_p1 = grid.index(x+1, y);
    //   atlas::idx_t idx_m1 = grid.index(x-1, y);

    //   if (halo_view(idx_p1)) { std::cout << "DBG " << idx << " " << idx_p1<<std::endl;}
    //   // // atlas::idx_t idx_p1 = idx+1;
    //   // // atlas::idx_t idx_m1 = idx-1;
    //   // if(halo_view(idx+1)) continue;

    //   // v_p1_view(idx, 0) = idx_p1 - idx; // v_view(idx_p1, 0);
    // }
    // geom_.functionSpace().haloExchange(fs[v]);
    // geom_.functionSpace().haloExchange(fs_p1[v]);
  }

  xx.validTime() += tstep_;
}

// -----------------------------------------------------------------------------

  void ModelSimpleWave::finalize(State & xx) const {
    fs_m1_.clear();
  }

  // -----------------------------------------------------------------------------
}  // namespace genericMarine
