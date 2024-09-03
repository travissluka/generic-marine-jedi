
/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "genericMarine/Traits.h"
#include "genericMarine/Model/ModelZonalAdvection.h"

#include "oops/interface/ModelBase.h"
#include "oops/util/FieldSetHelpers.h"
namespace genericMarine {

// -----------------------------------------------------------------------------

static oops::interface::ModelMaker<Traits, ModelZonalAdvection> modelAdvection_("ZonalAdvection");

// -----------------------------------------------------------------------------

double interp(const double lat, const std::vector<double> & lats,
              const std::vector<double> & vals ) {
  // find index i, and i-1 that are involved in the linear interp
  size_t i = 0;
  for (i = 0; i < lats.size(); ++i) {
    if (lat >= lats[i]) break;
  }

  // interpolation
  if (i == 0) {
    return vals[0];
  } else if (i == lats.size()) {
    return vals[vals.size()-1];
  } else {
    return (lat - lats[i])/(lats[i-1]-lats[i]) * (vals[i-1] - vals[i]) + vals[i];
  }
}

// -----------------------------------------------------------------------------

ModelZonalAdvection::ModelZonalAdvection(const Geometry & geom,
                                         const eckit::Configuration & conf)
  : ModelAdvectionBase(geom, oops::validateAndDeserialize<ModelAdvectionBase::Parameters>(conf))
{
  ModelZonalAdvection::Parameters params  = oops::validateAndDeserialize<ModelZonalAdvection::Parameters>(conf);

  // get advection parameters
  ModelZonalAdvection::Parameters::Advection advec_param =
    oops::validateAndDeserialize<ModelZonalAdvection::Parameters::Advection>(params.advection.value().config.value());
  const auto & speed_param = advec_param.speed.value();
  const std::vector<double> & lats = speed_param.latitude.value();
  const std::vector<double> & speed = speed_param.value.value();
  ASSERT(lats.size() == speed.size());
  ASSERT(lats.size() >= 1);
  for (size_t i = 1; i < lats.size(); i++) {
    ASSERT(lats[i] < lats[i-1]); }

  // get fields
  atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  auto cx = atlas::array::make_view<double, 1> (phaseSpeed_.field("cx"));
  auto cy = atlas::array::make_view<double, 1> (phaseSpeed_.field("cy"));
  auto coastdist = atlas::array::make_view<double, 2>(geom_.fields().field("distanceToCoast"));
  auto lonlat = atlas::array::make_view<double, 2>(fspace.lonlat());

  // set a horizontally varying u
  const auto & coast_damp_param = advec_param.coastDamp.value();
  const double coast_dist = coast_damp_param.distance.value();
  const double coast_damp = coast_damp_param.damping.value();
  for (atlas::idx_t idx = 0; idx < fspace.size(); idx++) {
    cx(idx) = interp(lonlat(idx, 1), lats, speed);
    if (coast_damp > 0.0 && coast_dist > 0.0) {
      double  d = 1.0 - ((1.0-coastdist(idx, 0)/coast_dist) * coast_damp);
      cx(idx) *= std::min(d, 1.0);
    }
  }

  // update the shear based diffusion
  updateDiffusionParams();

}

// -----------------------------------------------------------------------------

}  // namespace genericMarine
