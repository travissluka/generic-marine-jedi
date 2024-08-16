/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/Parameter.h"

// forward declarations
namespace genericMarine {
  class Geometry;
  class ModelAuxControl;
  class State;
  struct Traits;
}

namespace genericMarine {

//-----------------------------------------------------------------------------

  class BoundaryConditionParameters:public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(BoundaryConditionParameters, Parameters)
   public:
    // Dirichlet boundary conditions for incoming flow.
    // Value at boundary = a*f_x0 + b where f_x0 is the value of a neighboring valid grid point.
    // Outflow assumes Neumann conditions.
    oops::RequiredParameter<double> a{"a", this};
    oops::RequiredParameter<double> b{"b", this};
  };

  class ModelAdvectionBaseParameters:public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(ModelAdvectionBaseParameters, Parameters)
   public:
    oops::OptionalParameter<std::string> name{"name", this};
    oops::RequiredParameter<util::Duration> tstep{"tstep", this};
    oops::RequiredParameter<BoundaryConditionParameters> boundary{"boundary condition", this};
    oops::RequiredParameter<oops::Variables> vars{"variables", this};
    oops::Parameter<double> asselinFilter{"asselin filter", 0.2, this};
    oops::Parameter<double> diffusion{"diffusion", 0.0, this};
  };

//-----------------------------------------------------------------------------

  class ModelAdvectionBase:public oops::interface::ModelBase<Traits>,
              private util::ObjectCounter<ModelAdvectionBase>
  {
   public:
    typedef ModelAdvectionBaseParameters Parameters_;
    ModelAdvectionBase(const Geometry &, const ModelAdvectionBaseParameters &);
    virtual ~ModelAdvectionBase() = 0;

    // main model run methods
    void initialize(State &) const;
    void step(State &, const ModelAuxControl &) const;
    void finalize(State &) const;

    // accessors
    const util::Duration & timeResolution() const {return tstep_;}
    const oops::Variables & variables() const {return vars_;}

   protected:
    const Geometry & geom_;
    atlas::FieldSet phaseSpeed_;

   private:
    void advectionStep(const atlas::Field &, atlas::Field &) const;
    void diffusionStep(const atlas::Field &, atlas::Field &) const;

    void print(std::ostream &) const;
    util::Duration tstep_;
    const oops::Variables vars_;

    // Value at boundary = bc_a*f_x0 + bc_b where
    // f_x0 is the value of a neighboring valid grid point.
    const double bc_a_, bc_b_;
    const double asselin_;
    const double diffusion_;
    mutable atlas::FieldSet xx_tm1_;  // model state at previous time, for leapfrog scheme
  };

//-----------------------------------------------------------------------------
}  // namespace genericMarine
