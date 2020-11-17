/*
 * (C) Copyright 2019-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_STATE_STATE_H_
#define UMDSST_STATE_STATE_H_

#include <string>

#include "umdsst/Fields/Fields.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

// forward declarations
namespace umdsst {
  class Geometry;
  class Increment;
}

namespace eckit {
  class Configuration;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}

// ----------------------------------------------------------------------------

namespace umdsst {

  // State class
  class State : private util::ObjectCounter<State>,
                public umdsst::Fields {
   public:
    static const std::string classname() {return "umdsst::State";}

    // constructors, destructors
    State(const Geometry &, const eckit::Configuration &);
    State(const Geometry &, const oops::Variables &,
          const util::DateTime &);
    State(const Geometry &, const State &);
    State(const State &);
    ~State();

    // wrappers of methods that are fully implemented in Fields
    State & operator+=(const Increment &);
  };
}  // namespace umdsst

#endif  // UMDSST_STATE_STATE_H_
