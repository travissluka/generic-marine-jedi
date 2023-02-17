/*
 * (C) Copyright 2019-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "genericMarine/Covariance/Covariance.h"
#include "genericMarine/Geometry/Geometry.h"
#include "genericMarine/Increment/Increment.h"
#include "genericMarine/LinearVariableChange/LinearVariableChange.h"
#include "genericMarine/ModelAux/ModelAuxControl.h"
#include "genericMarine/ModelAux/ModelAuxCovariance.h"
#include "genericMarine/ModelAux/ModelAuxIncrement.h"
#include "genericMarine/State/State.h"
#include "genericMarine/VariableChange/VariableChange.h"

namespace genericMarine {

  struct Traits{
    static std::string name() {return "genericMarine";}
    static std::string nameCovar() {return "genericMarineCovar";}
    static std::string nameCovar4D() {return "genericMarineCovar";}

    // Interfaces that genericMarine has to implement
    // ---------------------------------------------------
    typedef genericMarine::Covariance           Covariance;
    typedef genericMarine::Geometry             Geometry;
    typedef genericMarine::Increment            Increment;
    typedef genericMarine::LinearVariableChange LinearVariableChange;
    typedef genericMarine::ModelAuxControl      ModelAuxControl;
    typedef genericMarine::ModelAuxCovariance   ModelAuxCovariance;
    typedef genericMarine::ModelAuxIncrement    ModelAuxIncrement;
    typedef genericMarine::State                State;
    typedef genericMarine::VariableChange       VariableChange;
  };
}  // namespace genericMarine

