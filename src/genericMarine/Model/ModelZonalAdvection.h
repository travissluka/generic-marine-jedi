
/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "genericMarine/Model/ModelAdvectionBase.h"

#include "oops/util/parameters/Parameters.h"

namespace genericMarine {

//-----------------------------------------------------------------------------

  class SpeedParameter:public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(SpeedParameter, oops::Parameters)
   public:
    oops::RequiredParameter<std::vector<double> > latitude{"latitude", this};
    oops::RequiredParameter<std::vector<double> > value{"value", this};
  };

  class ModelZonalAdvectionParameters:public ModelAdvectionBaseParameters {
    OOPS_CONCRETE_PARAMETERS(ModelZonalAdvectionParameters, ModelAdvectionBaseParameters)
   public:
    oops::RequiredParameter<SpeedParameter> speed{"speed", this};
    oops::RequiredParameter<double> coastDist{"coastal damping distance", this};
  };

//-----------------------------------------------------------------------------

  class ModelZonalAdvection:public ModelAdvectionBase,
    private util::ObjectCounter<ModelZonalAdvection> {
   public:
    typedef ModelZonalAdvectionParameters Parameters_;
    ModelZonalAdvection(const Geometry &, const ModelZonalAdvectionParameters &);
  };

//-----------------------------------------------------------------------------
}  // namespace genericMarine
