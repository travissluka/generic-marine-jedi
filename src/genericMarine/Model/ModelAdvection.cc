/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Traits.h"
#include "genericMarine/Model/ModelAdvection.h"
#include "genericMarine/Geometry/Geometry.h"

#include "oops/interface/ModelBase.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/missingValues.h"

#include "oops/util/abor1_cpp.h"

namespace genericMarine {

// -----------------------------------------------------------------------------

ModelAdvection::ModelAdvection(const Geometry & geom, const ModelAdvectionParameters & params)
: geom_(geom), tstep_(params.tstep), phaseSpeed_() {

  // create zero u/v fields
  atlas::Field cx = geom_.functionSpace().createField<double>(atlas::option::name("cx"));
  atlas::Field cy = geom_.functionSpace().createField<double>(atlas::option::name("cy"));
  auto cx_view = atlas::array::make_view<double, 1> (cx);
  auto cy_view = atlas::array::make_view<double, 1> (cy);
  cx_view.assign(0.0);
  cy_view.assign(0.0);
  phaseSpeed_.add(cx);
  phaseSpeed_.add(cy);
}

// -----------------------------------------------------------------------------

ModelAdvection::~ModelAdvection(){}

// -----------------------------------------------------------------------------

void ModelAdvection::print(std::ostream & os) const {
  os << "ModelAdvection::print not implemented";
}

// -----------------------------------------------------------------------------

void ModelAdvection::initialize(State & xx) const {
  xx_tm1_.clear();
}

// -----------------------------------------------------------------------------

// atlas::FieldSet ModelAdvection::setParams() const {
//   // This is a simple test situation, of 0 speed at the poles, increasing toward the equator
//   // with a modulation to keep near 0.0 at the coasts
//   atlas::FieldSet fset;

//   // other fields we'll need
//   atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
//   auto v_lonlat = atlas::array::make_view<double, 2>(fspace.lonlat());
//   auto v_halo = atlas::array::make_view<int, 1>(fspace.ghost());
//   auto v_coastdist = atlas::array::make_view<double, 2>(geom_.extraFields().field("distanceToCoast"));

//   // create zero u/v fields
//   atlas::Field cx = fspace.createField<double>(atlas::option::name("cx"));
//   atlas::Field cy = fspace.createField<double>(atlas::option::name("cy"));
//   auto cx_view = atlas::array::make_view<double, 1> (cx);
//   auto cy_view = atlas::array::make_view<double, 1> (cy);
//   cx_view.assign(0.0);
//   cy_view.assign(0.0);
//   fset.add(cx);
//   fset.add(cy);

//   // set a horizontally varying u
//   const double lat0 = 80.0;
//   const double lat1_val = -1.0;
//   const double coast_dist = 300e3;
//   for(atlas::idx_t idx = 0; idx < fspace.size(); idx++){
//     // based on latitude
//     double lat = v_lonlat(idx, 1);
//     if (abs(lat) > lat0) {
//       cx_view(idx) = 0.0;
//     } else {
//       cx_view(idx) = (1.0 - abs(lat)/lat0) * lat1_val;
//     }

//     // set to zero near coast
//     cx_view(idx) *= std::min(v_coastdist(idx,0) / coast_dist, 1.0);
//   }

//   // keep horizontally varying v to 0
//   return fset;
// }

// -----------------------------------------------------------------------------

void ModelAdvection::step(State & xx, const ModelAuxControl &) const {
  atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  double missing; missing = util::missingValue(missing);

  // get various data views we need
  auto dx = atlas::array::make_view<double, 2>(geom_.extraFields().field("dx"));
  auto dy = atlas::array::make_view<double, 2>(geom_.extraFields().field("dy"));
  auto cx = atlas::array::make_view<double, 1>(phaseSpeed_.field("cx"));
  auto cy = atlas::array::make_view<double, 1>(phaseSpeed_.field("cy"));

  // fieldsets at various times (past, present, future)
  atlas::FieldSet xx_tp1 = xx.fieldSet();
  atlas::FieldSet xx_t0 = util::copyFieldSet(xx_tp1);
  bool leapfrog_init = xx_tm1_.empty();
  if (leapfrog_init) xx_tm1_ = util::copyFieldSet(xx_t0);

  // temporary working fields, used later
  atlas::Field dfdx_field = fspace.createField<double>();
  atlas::Field dfdy_field = fspace.createField<double>();
  auto dfdx = atlas::array::make_view<double, 1>(dfdx_field);
  auto dfdy = atlas::array::make_view<double, 1>(dfdy_field);

  // for each variable in the model
  for (atlas::idx_t var = 0; var < xx_t0.size(); var++) {
    // various data views we need
    auto f_t0 = atlas::array::make_view<double,2>(xx_t0[var]);
    auto f_tp1 = atlas::array::make_view<double,2>(xx_tp1[var]);
    auto f_tm1 = atlas::array::make_view<double,2>(xx_tm1_[var]);

    // calculate the horizontal derivative
    dfdx.assign(0.0);
    dfdy.assign(0.0);
    for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
      for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
        atlas::idx_t idx = fspace.index(ii, jj);
        atlas::idx_t idx_xp1 = fspace.index(ii+1, jj);
        atlas::idx_t idx_xm1 = fspace.index(ii-1, jj);
        atlas::idx_t idx_yp1 = fspace.index(ii, jj+1);
        atlas::idx_t idx_ym1 = fspace.index(ii, jj-1);

        // skip land
        if (f_t0(idx,0) == missing) continue;

        // df/dx
        if (f_t0(idx_xp1, 0) == missing && f_t0(idx_xm1, 0) == missing) {
          // leave as 0, no neighbors
        } else if (f_t0(idx_xm1, 0) == missing) {
          // no left neighbor
          dfdx(idx) = (f_t0(idx_xp1, 0) - f_t0(idx, 0)) / dx(idx, 0);
        } else if (f_t0(idx_xp1, 0) == missing) {
          // no right neighbor
          dfdx(idx) = (f_t0(idx, 0) - f_t0(idx_xm1, 0)) / dx(idx, 0);
        } else {
          // 2 neighbors, centered difference
          dfdx(idx) = (f_t0(idx_xp1, 0) - f_t0(idx_xm1, 0)) / (2.0*dx(idx, 0));
        }

        // df/dy
        if (f_t0(idx_yp1, 0) == missing && f_t0(idx_ym1, 0) == missing) {
          // leave as 0, no neighbors
        } else if (f_t0(idx_ym1, 0) == missing) {
          // no neighbor below
          dfdy(idx) = (f_t0(idx_yp1, 0) - f_t0(idx, 0)) / dy(idx, 0);
        } else if (f_t0(idx_yp1, 0) == missing) {
          // no neighbor above
          dfdy(idx) = (f_t0(idx, 0) - f_t0(idx_ym1, 0)) / dy(idx, 0);
        } else {
          // 2 neighbors, centered difference
          dfdy(idx) = (f_t0(idx_yp1, 0) - f_t0(idx_ym1, 0)) / (2.0*dy(idx, 0));
        }
        // TODO double check these signs, it depends on how grid is loaded which way is up
        // this assumes positive lat is at start of grid (opposite of what you'd expect)
        dfdy(idx) *= -1.0;
      }
    }

    // time derivatives
    for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
      for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
        atlas::idx_t idx = fspace.index(ii, jj);

        // skip land
        if (f_t0(idx,0) == missing) continue;

        // calculate and apply add df/dt
        double dfdt = cx(idx) * dfdx(idx) + cy(idx) * dfdy(idx);
        if (leapfrog_init) {
          // euler forward (for very first timestep only)
          f_tp1(idx, 0) = f_t0(idx, 0) - tstep_.toSeconds()*dfdt ;
        } else {
          // leapfrog
          f_tp1(idx, 0) = f_tm1(idx, 0) - 2.0*tstep_.toSeconds()*dfdt;
        }

        // move t=0 to t-1
        f_tm1(idx, 0) = f_t0(idx, 0);
      }
    }

    // update halos
    fspace.haloExchange(xx_tp1[var]);
    fspace.haloExchange(xx_tm1_[var]);
  }

  xx.validTime() += tstep_;
}

// -----------------------------------------------------------------------------

  void ModelAdvection::finalize(State & xx) const {
    xx_tm1_.clear();
  }

  // -----------------------------------------------------------------------------
}  // namespace genericMarine
