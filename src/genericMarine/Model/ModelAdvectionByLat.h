
/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "genericMarine/Model/ModelAdvection.h"

namespace genericMarine {

//-----------------------------------------------------------------------------

  class ModelAdvectionByLatParameters:public ModelAdvectionParameters {
    OOPS_CONCRETE_PARAMETERS(ModelAdvectionByLatParameters, ModelAdvectionParameters)
   public:
    oops::RequiredParameter<std::vector<double> > latitude{"latitude", this};
    oops::RequiredParameter<std::vector<double> > phaseSpeed{"phase speed", this};
  };

//-----------------------------------------------------------------------------

  class ModelAdvectionByLat:public ModelAdvection {
   public:
    typedef ModelAdvectionByLatParameters Parameters_;
    ModelAdvectionByLat(const Geometry &, const ModelAdvectionByLatParameters &);
  };

//-----------------------------------------------------------------------------
}