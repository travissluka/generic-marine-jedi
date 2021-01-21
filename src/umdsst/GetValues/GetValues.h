/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_GETVALUES_GETVALUES_H_
#define UMDSST_GETVALUES_GETVALUES_H_

#include <memory>
#include <ostream>
#include <string>

#include "umdsst/GetValues/LocationsWrapper.h"
#include "umdsst/GetValues/GeoVaLsWrapper.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/Locations.h"

// forward declarations
namespace umdsst {
  class Geometry;
  class State;
}

namespace eckit {
  class Configuration;
}
namespace oops {
  class InterpolatorUnstructured;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}


// ----------------------------------------------------------------------------

namespace umdsst {

  // GetValues class: interpolate state to observation locations
  class GetValues : public util::Printable,
                    private util::ObjectCounter<GetValues> {
   public:
    static const std::string classname() {return "umdsst::GetValues";}

    // constructors, destructors
    GetValues(const Geometry &,
              const ufo::Locations & locs,
              const eckit::Configuration & config);
    virtual ~GetValues();

    // fills in geovals for all observations in the timeframe (t1, t2],
    void fillGeoVaLs(const State &,
                     const util::DateTime & t1,
                     const util::DateTime & t2,
                     ufo::GeoVaLs &) const;

   private:
    void print(std::ostream &) const;

    std::unique_ptr<oops::InterpolatorUnstructured> interpolator_;
    std::shared_ptr<const Geometry> geom_;
    LocationsWrapper locs_;
  };
}  // namespace umdsst

#endif  // UMDSST_GETVALUES_GETVALUES_H_
