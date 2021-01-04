/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "umdsst/GetValues/LinearGetValues.h"
#include "umdsst/State/State.h"
#include "umdsst/Geometry/Geometry.h"
#include "umdsst/Increment/Increment.h"
#include "umdsst/GetValues/GeoVaLsWrapper.h"
#include "umdsst/GetValues/GeoVaLsWrapperAD.h"

#include "eckit/config/Configuration.h"

#include "oops/generic/InterpolatorUnstructured.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace umdsst {

// ----------------------------------------------------------------------------

  LinearGetValues::LinearGetValues(const Geometry & geom,
                                   const ufo::Locations & locs)
    : geom_( new Geometry(geom)), locs_(locs) {
    interpolator_.reset( new oops::InterpolatorUnstructured(
                               eckit::LocalConfiguration(),
                               *geom_->atlasFunctionSpace(),
                               locs_.atlasFunctionSpace()) );
  }

// ----------------------------------------------------------------------------

  LinearGetValues::~LinearGetValues() { }

// ----------------------------------------------------------------------------

  void LinearGetValues::fillGeoVaLsAD(Increment & inc,
                                      const util::DateTime & t1,
                                      const util::DateTime & t2,
                                      const ufo::GeoVaLs & geovals) const {
    oops::Variables vars = geovals.getVars();
    for (size_t i = 0; i < vars.size(); i++) {
      // expect only sst for now
      if (vars[i] != "sea_surface_temperature")
        util::abor1_cpp("LinearGetValues::fillGeoVaLsAD,unkown state variable");

      atlas::Field fgvl = locs_.atlasFunctionSpace()->createField<double>(
                                                 atlas::option::levels(1));
      // copy from geovals to fin so it can be used in apply_ad;
      GeoVaLsWrapperAD(geovals, locs_.locs()).fill(t1, t2, fgvl);

      interpolator_->apply_ad(fgvl,
        inc.atlasFieldSet()->field("sea_surface_temperature"));
    }
  }

// ----------------------------------------------------------------------------

  void LinearGetValues::fillGeoVaLsTL(const Increment & inc,
                                      const util::DateTime & t1,
                                      const util::DateTime & t2,
                                      ufo::GeoVaLs & geovals) const {
    oops::Variables vars = geovals.getVars();
    for (size_t i = 0; i < vars.size(); i++) {
      // expect only sst for now
      if (vars[i] != "sea_surface_temperature")
        util::abor1_cpp("LinearGetValues::fillGeoVaLsTL,unkown state variable");

      atlas::Field fout = locs_.atlasFunctionSpace()->createField<double>(
                                                 atlas::option::levels(1));

      interpolator_->apply(
        inc.atlasFieldSet()->field("sea_surface_temperature"), fout);

      GeoVaLsWrapper(geovals, locs_.locs()).fill(t1, t2, fout);
    }
  }

// ----------------------------------------------------------------------------

  void LinearGetValues::setTrajectory(const State & state,
                                      const util::DateTime & t1,
                                      const util::DateTime & t2,
                                      ufo::GeoVaLs & geovals) {
    oops::Variables vars = geovals.getVars();
    for (size_t i = 0; i < vars.size(); i++) {
      // expect only sst for now
      if (vars[i] != "sea_surface_temperature")
        util::abor1_cpp("LinearGetValues::setTrajectory,unkown state variable");

      atlas::Field fout = locs_.atlasFunctionSpace()->createField<double>(
                                                 atlas::option::levels(1));

      interpolator_->apply(
        state.atlasFieldSet()->field("sea_surface_temperature"), fout);

      GeoVaLsWrapper(geovals, locs_.locs()).fill(t1, t2, fout);
    }
  }

// ----------------------------------------------------------------------------

  void LinearGetValues::print(std::ostream & os) const {
    os << "umdsst::LinearGetValues" << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace umdsst
