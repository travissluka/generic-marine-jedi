/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Traits.h"
#include "genericMarine/Model/ModelAdvectionBase.h"
#include "genericMarine/Geometry/Geometry.h"

#include "oops/interface/ModelBase.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/missingValues.h"

#include "oops/util/abor1_cpp.h"

namespace genericMarine {

// -----------------------------------------------------------------------------

ModelAdvectionBase::ModelAdvectionBase(const Geometry & geom, const ModelAdvectionBaseParameters & params)
 : geom_(geom), tstep_(params.tstep), phaseSpeed_(), 
   bc_a_(params.boundary.value().a),
   bc_b_(params.boundary.value().b) {

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

ModelAdvectionBase::~ModelAdvectionBase(){}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::print(std::ostream & os) const {
  os << "ModelAdvectionBase::print not implemented";
}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::initialize(State & xx) const {
  xx_tm1_.clear();
}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::step(State & xx, const ModelAuxControl &) const {
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

        // note on boundary conditions for df/dx and df/dy:
        // For outgoing flow, Neumann conditions are used (derivative is specified) so that the second derivative is 0
        // For incoming flow, boundary value = bc_a_*f_x0 + bc_b_ where f_x0 is a valid neighboring grid point
        double f_x0  = f_t0(idx, 0);
        double inflow_bc = bc_a_*f_x0 + bc_b_;

        // df/dx
        double f_xm1 = f_t0(idx_xm1, 0);        
        double f_xp1 = f_t0(idx_xp1, 0);
        if (f_xm1 == missing){
          // set boundary condition for missing left neighbor
          f_xm1 = inflow_bc;
          if (cx(idx) < 0.0 && f_xp1 != missing) f_xm1 = 2.0*f_x0 - f_xp1;
        } 
        if (f_xp1 == missing){
        // set boundary condition for missing right neighbor
         f_xp1 = inflow_bc;
         if (cx(idx) > 0.0 && f_xm1 != missing) f_xp1 =  2.0*f_x0 - f_xm1;
        }
        dfdx(idx) = (f_xp1 - f_xm1) / (2.0*dx(idx, 0));

        // df/dy
        double f_ym1 = f_t0(idx_ym1, 0);
        double f_yp1 = f_t0(idx_yp1, 0);
        if (f_ym1 == missing){
          // set boundary condition for missing neighbor below
          f_ym1 = inflow_bc;
          if (cy(idx) < 0.0 && f_yp1 != missing) f_ym1 = 2.0*f_x0 - f_yp1;
        }
        if (f_yp1 == missing) {
          // set boundary condition for missing neighbor above
          f_yp1 = inflow_bc;
          if (cy(idx) > 0.0 && f_ym1 != missing) f_yp1 = 2.0*f_x0 - f_ym1;
        }
        dfdy(idx) = (f_yp1 - f_ym1) / (2.0*dy(idx, 0));     
        // TODO double check the signs for cy, it depends on how grid is loaded which way is up
        // this assumes positive lat is at start of grid (opposite of what you'd expect)
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

  void ModelAdvectionBase::finalize(State & xx) const {
    xx_tm1_.clear();
  }

  // -----------------------------------------------------------------------------
}  // namespace genericMarine
