/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"

// forward declarations
namespace genericMarine {
  class Geometry;
  class ModelAuxControl;
  class State;
  struct Traits;
}

namespace genericMarine {

//-----------------------------------------------------------------------------

  class ModelSimpleWaveParameters:public oops::ModelParametersBase {
    OOPS_CONCRETE_PARAMETERS(ModelSimpleWaveParameters, ModelParametersBase)
   public:
    oops::RequiredParameter<util::Duration> tstep{"tstep", this};
  };

//-----------------------------------------------------------------------------

  class ModelSimpleWave:public oops::interface::ModelBase<Traits>,
              private util::ObjectCounter<ModelSimpleWave>
  {
   public:
    typedef ModelSimpleWaveParameters Parameters_;
    ModelSimpleWave(const Geometry &, const ModelSimpleWaveParameters &);

    // main model run methods
    void initialize(State &) const;
    void step(State &, const ModelAuxControl &) const;
    void finalize(State &) const;

    // accessors
    const util::Duration & timeResolution() const {return tstep_;}
    const oops::Variables & variables() const {return vars_;}

   private:
    void print(std::ostream &) const;
    util::Duration tstep_;
    const oops::Variables vars_;
    const Geometry & geom_;

    mutable atlas::FieldSet fs_m1_;
    atlas::Field dx_;
  };

//-----------------------------------------------------------------------------
}  // namespace genericMarine
