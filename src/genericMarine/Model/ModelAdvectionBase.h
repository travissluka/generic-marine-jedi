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
#include "oops/util/parameters/ConfigurationParameter.h"

// forward declarations
namespace genericMarine {
  class Geometry;
  class ModelAuxControl;
  class State;
  struct Traits;
}

namespace genericMarine {

//-----------------------------------------------------------------------------

  class ModelAdvectionBase:public oops::interface::ModelBase<Traits>,
              private util::ObjectCounter<ModelAdvectionBase>
  {
   public:
    // -----------------------------------------------------------------------------
    class Parameters:public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Parameters, oops::Parameters)

     public:
       // -----------------------------------------------------------------------------
      class BoundaryCondition:public oops::Parameters {
        OOPS_CONCRETE_PARAMETERS(BoundaryCondition, oops::Parameters)
       public:
        // Dirichlet boundary conditions for incoming flow.
        // Value at boundary = a*f_x0 + b where f_x0 is the value of a neighboring valid grid point.
        // Outflow assumes Neumann conditions.
        oops::RequiredParameter<double> a{"a", this};
        oops::RequiredParameter<double> b{"b", this};
      };
      // -----------------------------------------------------------------------------
      class AdvectionBase:public oops::Parameters {
        OOPS_ABSTRACT_PARAMETERS(AdvectionBase, oops::Parameters)
       public:
        oops::RequiredParameter<BoundaryCondition> boundary{"boundary condition", this};
        oops::Parameter<double> asselinFilter{"asselin filter", 0.2, this};
      };
      class Advection:public AdvectionBase {
        OOPS_CONCRETE_PARAMETERS(Advection, AdvectionBase)
       public:
        oops::ConfigurationParameter config{this};
      };

      // -----------------------------------------------------------------------------
      class Diffusion:public oops::Parameters {
        OOPS_CONCRETE_PARAMETERS(Diffusion, oops::Parameters)
       public:
        oops::Parameter<int> smootherIterations{"coefficient smoothing", 1, this};
        oops::Parameter<double> kh{"Kh", 0.0, this};
        oops::Parameter<double> ah{"Ah", 0.0, this};
        oops::Parameter<double> kh_smag{"Kh_smag scale", 0.0, this};
        oops::Parameter<double> kh_smag_max{"Kh_smag max", 0.0, this};
      };
      // -----------------------------------------------------------------------------
      oops::OptionalParameter<std::string> name{"name", this};
      oops::RequiredParameter<util::Duration> tstep{"tstep", this};
      oops::RequiredParameter<oops::Variables> vars{"variables", this};
      oops::RequiredParameter<Advection> advection{"advection", this};
      oops::Parameter<Diffusion> diffusion{"diffusion", {}, this};
    };
    // -----------------------------------------------------------------------------

    typedef Parameters Parameters_;
    ModelAdvectionBase(const Geometry &, const Parameters &);
    virtual ~ModelAdvectionBase() = 0;

    // main model run methods
    void initialize(State &) const;
    void step(State &, const ModelAuxControl &) const;
    void finalize(State &) const;

    // accessors
    const util::Duration & timeResolution() const {return tstep_;}
    const oops::Variables & variables() const {return vars_;}

   protected:
    void updateDiffusionParams();

    const Geometry & geom_;
    atlas::FieldSet phaseSpeed_;
    Parameters::Diffusion diffusionParams_;
    atlas::FieldSet diffusionCoeffs_;

   private:
    void advectionStep(const atlas::Field &, atlas::Field &) const;
    void diffusionStep(const atlas::Field &, atlas::Field &, double) const;

    void print(std::ostream &) const;
    util::Duration tstep_;
    const oops::Variables vars_;

    // Value at boundary = bc_a*f_x0 + bc_b where
    // f_x0 is the value of a neighboring valid grid point.
    const double bc_a_, bc_b_;
    const double asselin_;
    mutable atlas::FieldSet xx_tm1_;  // model state at previous time, for leapfrog scheme
  };

//-----------------------------------------------------------------------------
}  // namespace genericMarine
