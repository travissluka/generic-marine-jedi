/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_TRAITS_H_
#define UMDSST_TRAITS_H_

#include <string>

#include "umdsst/Covariance/Covariance.h"
#include "umdsst/Geometry/Geometry.h"
#include "umdsst/GetValues/GetValues.h"
#include "umdsst/GetValues/LinearGetValues.h"
#include "umdsst/Increment/Increment.h"
#include "umdsst/ModelAux/ModelAuxControl.h"
#include "umdsst/ModelAux/ModelAuxCovariance.h"
#include "umdsst/ModelAux/ModelAuxIncrement.h"
#include "umdsst/State/State.h"

namespace umdsst {

  struct Traits{
    static std::string name() {return "umdsst";}
    static std::string nameCovar() {return "umdsstCovar";}
    static std::string nameCovar4D() {return "umdsstCovar";}

    // Interfaces that umdsst has to implement
    // ---------------------------------------------------
    typedef umdsst::Covariance          Covariance;
    typedef umdsst::Geometry            Geometry;
    typedef umdsst::GetValues           GetValues;
    typedef umdsst::Increment           Increment;
    typedef umdsst::LinearGetValues     LinearGetValues;
    typedef umdsst::ModelAuxControl     ModelAuxControl;
    typedef umdsst::ModelAuxCovariance  ModelAuxCovariance;
    typedef umdsst::ModelAuxIncrement   ModelAuxIncrement;
    typedef umdsst::State               State;
  };
}  // namespace umdsst

#endif  // UMDSST_TRAITS_H_
