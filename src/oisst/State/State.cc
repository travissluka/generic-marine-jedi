/*
 * (C) Copyright 2019-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "oisst/Geometry/Geometry.h"
#include "oisst/Increment/Increment.h"
#include "oisst/State/State.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oisst {

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const eckit::Configuration & conf)
    : Fields(geom, oops::Variables(conf, "state variables"),
             util::DateTime(conf.getString("date"))) {
    // check and read in data if it has filename
    if (conf.has("filename"))
      read(conf);
  }

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & time)
    : Fields(geom, vars, time) {}

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const State & other)
    : State(other) {
    // Change state resolution, normally used for interpolation.
    // For now, just copy without interpolation
  }

// ----------------------------------------------------------------------------

  State::State(const State & other)
    : Fields(other) {}

// ----------------------------------------------------------------------------

  State::~State() {}

// ----------------------------------------------------------------------------

  State & State::operator+=(const Increment & dx) {
    Fields::operator+=(dx);
    return *this;
  }

// ----------------------------------------------------------------------------

}  // namespace oisst
