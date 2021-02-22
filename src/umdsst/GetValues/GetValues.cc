/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "umdsst/Geometry/Geometry.h"
#include "umdsst/GetValues/GetValues.h"
#include "umdsst/GetValues/GeoVaLsWrapper.h"
#include "umdsst/State/State.h"

#include "eckit/config/Configuration.h"

#include "oops/generic/InterpolatorUnstructured.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace umdsst {

// ----------------------------------------------------------------------------

  GetValues::GetValues(const Geometry & geom,
                       const ufo::Locations & locs,
                       const eckit::Configuration & config)
    : geom_(new Geometry(geom)), locs_(locs) {
    interpolator_.reset( new oops::InterpolatorUnstructured(
                               eckit::LocalConfiguration(),
                               *geom_->atlasFunctionSpace(),
                               locs_.atlasFunctionSpace()));
  }

// -----------------------------------------------------------------------------

  GetValues::~GetValues() { }

// ----------------------------------------------------------------------------

  void GetValues::fillGeoVaLs(const State & state,
                              const util::DateTime & t1,
                              const util::DateTime & t2,
                              ufo::GeoVaLs & geovals) const {
    oops::Variables vars = geovals.getVars();
    std::vector<atlas::Field> fields(vars.size());

    for (size_t i = 0; i < vars.size(); i++) {
      fields[i] = locs_.atlasFunctionSpace()->createField<double>(
                                                 atlas::option::levels(1));

      if (vars[i] == "sea_surface_temperature") {
        interpolator_->apply(
          state.atlasFieldSet()->field("sea_surface_temperature"), fields[i]);
      } else if (vars[i] == "sea_area_fraction") {
        interpolator_->apply(
          geom_->atlasFieldSet()->field("ggmask"), fields[i]);  // not gmask
      }
    }

    GeoVaLsWrapper(geovals, locs_.locs()).fill(t1, t2, fields);
  }

// ----------------------------------------------------------------------------

  void GetValues::print(std::ostream & os) const {
    os << "umdsst::GetValues" << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace umdsst
