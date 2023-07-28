
/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Traits.h"
#include "genericMarine/Model/ModelAdvectionByLat.h"

#include "oops/interface/ModelBase.h"

namespace genericMarine {

// -----------------------------------------------------------------------------

static oops::interface::ModelMaker<Traits, ModelAdvectionByLat> modelAdvection_("AdvectionByLat");

// -----------------------------------------------------------------------------

double interp(const double lat, const std::vector<double> & lats, const std::vector<double> & vals ) {
  size_t i = 0;
  for( i = 0; i < lats.size(); ++i) if (lat >= lats[i]) break;
  if (i == 0) {
    return vals[0];
  } else if (i == lats.size()) {
    return vals[vals.size()-1];
  } else {
    double a = (lat - lats[i])/(lats[i-1]-lats[i]);
    return a * (vals[i-1] - vals[i]) + vals[i];
  }
}

// -----------------------------------------------------------------------------

ModelAdvectionByLat::ModelAdvectionByLat(const Geometry & geom, const ModelAdvectionByLatParameters & params)
 : ModelAdvection(geom, params)
{
  // make sure input parameters are correct
  const std::vector<double> & lats = params.latitude.value();
  const std::vector<double> & vals = params.phaseSpeed.value();
  ASSERT (lats.size() == vals.size());
  ASSERT (lats.size() >= 1);
  for (size_t i = 1; i < lats.size(); i++) {
    ASSERT(lats[i] < lats[i-1]);
  }

  // get fields
  atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  auto cx = atlas::array::make_view<double, 1> (phaseSpeed_.field("cx"));
  auto cy = atlas::array::make_view<double, 1> (phaseSpeed_.field("cy"));
  auto coastdist = atlas::array::make_view<double, 2>(geom_.extraFields().field("distanceToCoast"));
  auto lonlat = atlas::array::make_view<double, 2>(fspace.lonlat());

  // set a horizontally varying u
  const double coast_dist = 300e3;
  for(atlas::idx_t idx = 0; idx < fspace.size(); idx++){
    // based on latitude
    cx(idx) = interp(lonlat(idx, 1), lats, vals);

    // set to zero near coast
    cx(idx) *= std::min(coastdist(idx,0) / coast_dist, 1.0);
  }
}

// -----------------------------------------------------------------------------

}