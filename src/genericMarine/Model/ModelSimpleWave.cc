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

namespace genericMarine {

  static oops::interface::ModelMaker<Traits, ModelSimpleWave> modelSimpleWave_("SimpleWave");

// -----------------------------------------------------------------------------

  ModelSimpleWave::ModelSimpleWave(const Geometry & geom, const ModelSimpleWaveParameters & params)
  : geom_(geom), tstep_(params.tstep) {
  }

// -----------------------------------------------------------------------------

  void ModelSimpleWave::print(std::ostream & os) const {
    os << "ModelSimpleWave::print not implemented";
  }

// -----------------------------------------------------------------------------

  void ModelSimpleWave::initialize(State & xx) const {
  }

// -----------------------------------------------------------------------------

  void ModelSimpleWave::step(State & xx, const ModelAuxControl &) const {
    xx.validTime() += tstep_;
  }

// -----------------------------------------------------------------------------

  void ModelSimpleWave::finalize(State & xx) const {
  }

  // -----------------------------------------------------------------------------
}  // namespace genericMarine
