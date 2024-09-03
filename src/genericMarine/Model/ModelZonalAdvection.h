
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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/NumericConstraints.h"

namespace eckit {
  class Configuration;
}

//-----------------------------------------------------------------------------

namespace genericMarine {

//-----------------------------------------------------------------------------
class ModelZonalAdvection:public ModelAdvectionBase,
  private util::ObjectCounter<ModelZonalAdvection> {
 public:
  // --------------------------------------------------------------------------
  class Parameters:public ModelAdvectionBase::Parameters {
    OOPS_CONCRETE_PARAMETERS(Parameters, ModelAdvectionBase::Parameters)
   public:
    // ------------------------------------------------------------------------
    class Speed:public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Speed, oops::Parameters)
     public:
      oops::RequiredParameter<std::vector<double> > latitude{"latitude", this};
      oops::RequiredParameter<std::vector<double> > value{"value", this};
    };
    // ------------------------------------------------------------------------
    class CoastalDamping:public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(CoastalDamping, oops::Parameters)
     public:
      oops::Parameter<double> distance{"distance", 0.0, this};
      oops::Parameter<double> damping{"amount", 1.0, this};
    };
    // ------------------------------------------------------------------------
    class Advection:public ModelAdvectionBase::Parameters::AdvectionBase {
      OOPS_CONCRETE_PARAMETERS(Advection, ModelAdvectionBase::Parameters::AdvectionBase)
     public:
      oops::RequiredParameter<Speed> speed{"speed", this};
      oops::Parameter<CoastalDamping> coastDamp{"coastal damping", {}, this};
    };

    // --------------------------------------------------------------------------
    // oops::RequiredParameter<Advection> advection{"advection", this};
  };

  // --------------------------------------------------------------------------

  typedef Parameters Parameters_;
  ModelZonalAdvection(const Geometry &, const eckit::Configuration &);
};

//-----------------------------------------------------------------------------
}  // namespace genericMarine
