/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "genericMarine/Traits.h"
#include "genericMarine/Model/ModelAdvectionBase.h"
#include "genericMarine/Geometry/Geometry.h"

#include "oops/interface/ModelBase.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/missingValues.h"

#include "oops/util/abor1_cpp.h"

namespace genericMarine {

// -----------------------------------------------------------------------------

ModelAdvectionBase::ModelAdvectionBase(const Geometry & geom,
                                       const Parameters & params)
  : geom_(geom), tstep_(params.tstep), phaseSpeed_(), vars_(params.vars),
    bc_a_(params.advection.value().boundary.value().a.value()),
    bc_b_(params.advection.value().boundary.value().b.value()),
    asselin_(params.advection.value().asselinFilter.value()),
    diffusionParams_(params.diffusion.value())
{
  // create zero u/v fields
  atlas::Field cx = geom_.functionSpace().createField<double>(atlas::option::name("cx"));
  atlas::Field cy = geom_.functionSpace().createField<double>(atlas::option::name("cy"));
  auto cx_view = atlas::array::make_view<double, 1> (cx);
  auto cy_view = atlas::array::make_view<double, 1> (cy);
  cx_view.assign(0.0);
  cy_view.assign(0.0);
  phaseSpeed_.add(cx);
  phaseSpeed_.add(cy);

  // initialize diffusion parameters
  atlas::Field kh = geom_.functionSpace().createField<double>(atlas::option::name("kh") | atlas::option::levels(1));
  atlas::Field ah = geom_.functionSpace().createField<double>(atlas::option::name("ah") | atlas::option::levels(1));
  diffusionCoeffs_.add(kh);
  diffusionCoeffs_.add(ah);
  updateDiffusionParams();
}

// -----------------------------------------------------------------------------

ModelAdvectionBase::~ModelAdvectionBase() {}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::print(std::ostream & os) const {
  os << "ModelAdvectionBase::print not implemented";
}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::initialize(State & xx) const {
  xx_tm1_.clear();
}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::updateDiffusionParams() {
  const atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  const double missing = util::missingValue<double>();
  const double maxDt = 2.0*tstep_.toSeconds();

  // get various views we'll need
  phaseSpeed_.haloExchange();
  const auto & cx = atlas::array::make_view<double, 1>(phaseSpeed_.field("cx"));
  const auto & cy = atlas::array::make_view<double, 1>(phaseSpeed_.field("cy"));
  const auto & dx = atlas::array::make_view<double, 2>(geom_.fields().field("dx"));
  const auto & dy = atlas::array::make_view<double, 2>(geom_.fields().field("dy"));
  const auto & mask = atlas::array::make_view<double, 2>(geom_.fields().field("mask"));
  auto & f_kh = diffusionCoeffs_.field("kh");
  auto & f_ah = diffusionCoeffs_.field("ah");
  auto kh = atlas::array::make_view<double, 2>(f_kh);
  auto ah = atlas::array::make_view<double, 2>(f_ah);

  // set to 0
  kh.assign(0.0);
  ah.assign(0.0);

  // calculate initial Kh, Ah
  for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
    for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
      const atlas::idx_t idx = fspace.index(ii, jj);
      const atlas::idx_t idx_xp1 = fspace.index(ii+1, jj);
      const atlas::idx_t idx_xm1 = fspace.index(ii-1, jj);
      const atlas::idx_t idx_yp1 = fspace.index(ii, jj+1);
      const atlas::idx_t idx_ym1 = fspace.index(ii, jj-1);

      // skip land
      if (mask(idx, 0) == 0.0 ) continue;

      // base diffusion
      kh(idx, 0) = diffusionParams_.kh;
      ah(idx, 0) = diffusionParams_.ah;

      // Smagorinsky horizontal diffusion
      // calculate shear (du/dy + dv/dx) and tension (du/dx - dv/dy)
      double shear = (cx(idx_yp1) - cx(idx_ym1)) / (2.0*dy(idx, 0)) +
                     (cy(idx_xp1) - cy(idx_xm1)) / (2.0*dx(idx, 0));
      double tension = (cx(idx_xp1) - cx(idx_xm1)) / (2.0*dx(idx, 0)) -
                       (cy(idx_yp1) - cy(idx_ym1)) / (2.0*dy(idx, 0));
      double kh_smag =  diffusionParams_.kh_smag * diffusionParams_.kh_smag * sqrt(tension*tension + shear*shear);

      // shear based diffusion
      kh(idx, 0) += std::min(kh_smag, 1.0 * diffusionParams_.kh_smag_max);
    }
  }

  // lambda function to smooth field
  auto smooth = [&](atlas::Field & field) {
    field.haloExchange();

    atlas::Field field_tmp = field.clone();
    const auto & field_tmp_view = atlas::array::make_view<double, 2>(field_tmp);
    auto field_view = atlas::array::make_view<double, 2>(field);

    for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
      for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
        const atlas::idx_t idx = fspace.index(ii, jj);
        const atlas::idx_t idx_xp1 = fspace.index(ii+1, jj);
        const atlas::idx_t idx_xm1 = fspace.index(ii-1, jj);
        const atlas::idx_t idx_yp1 = fspace.index(ii, jj+1);
        const atlas::idx_t idx_ym1 = fspace.index(ii, jj-1);

        // skip land
        if (mask(idx, 0) == 0.0 ) continue;

        int count = 1;
        double val = field_tmp_view(idx, 0);
        if (mask(idx_xp1, 0) > 0.0) { count++; val += field_tmp_view(idx_xp1, 0); }
        if (mask(idx_xm1, 0) > 0.0) { count++; val += field_tmp_view(idx_xm1, 0); }
        if (mask(idx_yp1, 0) > 0.0) { count++; val += field_tmp_view(idx_yp1, 0); }
        if (mask(idx_ym1, 0) > 0.0) { count++; val += field_tmp_view(idx_ym1, 0); }

        field_view(idx, 0) = val / count;
      }
    }
    field.set_dirty();
  };

  // smooth Kh, Ah
  const auto iterations = diffusionParams_.smootherIterations;
  for (int i = 0; i < iterations; i++) {
    smooth(f_kh);
    smooth(f_ah);
  }

  // clamp Kh, Ah
  for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
    for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
      const atlas::idx_t idx = fspace.index(ii, jj);

      // skip land
      if (mask(idx, 0) == 0.0 ) continue;

      // limit diffusion values where it would violate CFL condition
      auto minGrd = std::min(dx(idx, 0), dy(idx, 0));
      kh(idx, 0) = std::min(kh(idx, 0), 0.5*minGrd*minGrd / maxDt);
      ah(idx, 0) = std::min(ah(idx, 0), (1.0/16.0)*minGrd*minGrd*minGrd*minGrd / maxDt);
    }
  }

  diffusionCoeffs_.set_dirty();
}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::advectionStep(const atlas::Field & f, atlas::Field & tendency) const {
  const atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  const double missing = util::missingValue<double>();

  // get various views we'll need
  const auto & f_t0 = atlas::array::make_view<double, 2>(f);
  auto dfdt = atlas::array::make_view<double, 2>(tendency);
  const auto & dx = atlas::array::make_view<double, 2>(geom_.fields().field("dx"));
  const auto & dy = atlas::array::make_view<double, 2>(geom_.fields().field("dy"));
  const auto & cx = atlas::array::make_view<double, 1>(phaseSpeed_.field("cx"));
  const auto & cy = atlas::array::make_view<double, 1>(phaseSpeed_.field("cy"));
  phaseSpeed_.haloExchange();

  // calculate horizontal derivatives
  for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
    for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
      const atlas::idx_t idx = fspace.index(ii, jj);
      const atlas::idx_t idx_xp1 = fspace.index(ii+1, jj);
      const atlas::idx_t idx_xm1 = fspace.index(ii-1, jj);
      const atlas::idx_t idx_yp1 = fspace.index(ii, jj+1);
      const atlas::idx_t idx_ym1 = fspace.index(ii, jj-1);

      // skip land
      if (f_t0(idx, 0) == missing) continue;

      // note on boundary conditions for df/dx and df/dy:
      // For outgoing flow, Neumann conditions are used (derivative is specified) so that
      // the second derivative is 0.
      // For incoming flow, boundary value = bc_a_*f_x0 + bc_b_ where f_x0 is a valid
      // neighboring grid point
      const double f_x0  = f_t0(idx, 0);
      const double inflow_bc = bc_a_*f_x0 + bc_b_;

      // df/dx * dx/dt
      double f_xm1 = f_t0(idx_xm1, 0);
      double f_xp1 = f_t0(idx_xp1, 0);
      if (f_xm1 == missing) {
        // set boundary condition for missing left neighbor
        f_xm1 = inflow_bc;
        if (cx(idx) < 0.0 && f_xp1 != missing) f_xm1 = 2.0*f_x0 - f_xp1;
      }
      if (f_xp1 == missing) {
      // set boundary condition for missing right neighbor
       f_xp1 = inflow_bc;
       if (cx(idx) > 0.0 && f_xm1 != missing) f_xp1 =  2.0*f_x0 - f_xm1;
      }
      dfdt(idx, 0) -= cx(idx)* (f_xp1 - f_xm1) / (2.0*dx(idx, 0));

      // df/dy * dy/dt
      double f_ym1 = f_t0(idx_ym1, 0);
      double f_yp1 = f_t0(idx_yp1, 0);
      if (f_ym1 == missing) {
        // set boundary condition for missing neighbor below
        f_ym1 = inflow_bc;
        if (cy(idx) < 0.0 && f_yp1 != missing) f_ym1 = 2.0*f_x0 - f_yp1;
      }
      if (f_yp1 == missing) {
        // set boundary condition for missing neighbor above
        f_yp1 = inflow_bc;
        if (cy(idx) > 0.0 && f_ym1 != missing) f_yp1 = 2.0*f_x0 - f_ym1;
      }
      dfdt(idx, 0) -= cy(idx) * (f_yp1 - f_ym1) / (2.0*dy(idx, 0));
      // TODO(travis) double check the signs for cy, it depends on how grid is loaded which way
      // is up this assumes positive lat is at start of grid (opposite of what you'd expect)
    }
  }
  tendency.set_dirty();
}

const atlas::Field laplacian(const atlas::Field & f, const Geometry & geom) {
  const atlas::functionspace::StructuredColumns fspace(geom.functionSpace());
  const double missing = util::missingValue<double>();

  f.haloExchange();

  // create field and set to 0
  atlas::Field lap = f.clone();
  auto lap_view = atlas::array::make_view<double, 2>(lap);
  lap_view.assign(0.0);

  // get various views we'll need
  const auto & f_t0 = atlas::array::make_view<double, 2>(f);
  const auto & dx = atlas::array::make_view<double, 2>(geom.fields().field("dx"));
  const auto & dy = atlas::array::make_view<double, 2>(geom.fields().field("dy"));


  for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
    for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
      const atlas::idx_t idx = fspace.index(ii, jj);
      const atlas::idx_t idx_xp1 = fspace.index(ii+1, jj);
      const atlas::idx_t idx_xm1 = fspace.index(ii-1, jj);
      const atlas::idx_t idx_yp1 = fspace.index(ii, jj+1);
      const atlas::idx_t idx_ym1 = fspace.index(ii, jj-1);

      // skip land
      if (f_t0(idx, 0) == missing) continue;

      // derivatves
      double x = 0.0, y = 0.0;
      if (f_t0(idx_xp1, 0) != missing) x += f_t0(idx_xp1, 0) - f_t0(idx, 0);
      if (f_t0(idx_xm1, 0) != missing) x += f_t0(idx_xm1, 0) - f_t0(idx, 0);
      if (f_t0(idx_yp1, 0) != missing) y += f_t0(idx_yp1, 0) - f_t0(idx, 0);
      if (f_t0(idx_ym1, 0) != missing) y += f_t0(idx_ym1, 0) - f_t0(idx, 0);
      lap_view(idx, 0) = (x / (dx(idx, 0)*dx(idx, 0)) + y / (dy(idx, 0)*dy(idx, 0)));
    }
  }
  lap.set_dirty();
  return lap;
}

void ModelAdvectionBase::diffusionStep(const atlas::Field & f, atlas::Field & tendency,
     double max_dt) const
{
  const atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  const double missing = util::missingValue<double>();

  // get various views we'll need
  const auto & f_t0 = atlas::array::make_view<double, 2>(f);
  auto dfdt = atlas::array::make_view<double, 2>(tendency);
  const auto & dx = atlas::array::make_view<double, 2>(geom_.fields().field("dx"));
  const auto & dy = atlas::array::make_view<double, 2>(geom_.fields().field("dy"));
  const auto & kh = atlas::array::make_view<double, 2>(diffusionCoeffs_.field("kh"));
  const auto & ah = atlas::array::make_view<double, 2>(diffusionCoeffs_.field("ah"));

  // make sure halos are up to date
  f.haloExchange();
  diffusionCoeffs_.haloExchange();

  // calculate laplacians
  const auto lap = laplacian(f, geom_);
  const auto lap2 = laplacian(lap, geom_);
  const auto & lap_view = atlas::array::make_view<double, 2>(lap);
  const auto & lap2_view = atlas::array::make_view<double, 2>(lap2);

  for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
    for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
      const atlas::idx_t idx = fspace.index(ii, jj);

      // skip land
      if (f_t0(idx, 0) == missing) continue;

      // laplcian diffusion
      dfdt(idx, 0) += kh(idx, 0) * lap_view(idx, 0);

      // biharmonic diffusion
      dfdt(idx, 0) -= ah(idx, 0) * lap2_view(idx, 0);
    }
  }
  tendency.set_dirty();
}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::step(State & xx, const ModelAuxControl &) const {
  const atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());
  const double missing = util::missingValue<double>();
  const bool leapfrogInit = xx_tm1_.empty();
  if (leapfrogInit) xx_tm1_ = util::copyFieldSet(xx.fieldSet());

  // existing fields we'll need
  auto & f_t0 = xx.fieldSet()[0];  // TODO(travis) why is it missing variable names??
  auto v_t0 = atlas::array::make_view<double, 2>(f_t0);
  auto & f_tm1 = xx_tm1_[0];
  auto v_tm1 = atlas::array::make_view<double, 2>(f_tm1);

  // make sure halos are up to date
  f_t0.haloExchange();
  f_tm1.haloExchange();

  // initialize zero tendency (df/dt)
  atlas::Field dfdt = f_t0.clone();
  auto v_dfdt = atlas::array::make_view<double, 2>(dfdt);
  v_dfdt.assign(0.0);

  // calculate advection and diffusion tendencies
  advectionStep(f_t0, dfdt);
  diffusionStep(leapfrogInit ? f_t0 : xx_tm1_[0] , dfdt, 2.0*tstep_.toSeconds());

  // timestep with the tendencies
  const double dt = tstep_.toSeconds() * (leapfrogInit ? 1.0 : 2.0);
  const auto & f_prev = leapfrogInit ? f_t0 : f_tm1;
  const auto v_prev = atlas::array::make_view<double, 2>(f_prev);
  for (atlas::idx_t jj = fspace.j_begin(); jj < fspace.j_end(); jj++) {
    for (atlas::idx_t ii = fspace.i_begin(jj); ii < fspace.i_end(jj); ii++) {
      const atlas::idx_t idx = fspace.index(ii, jj);

      if (v_prev(idx, 0) == missing) continue;

      // apply timestep
      double v_tp1 = v_prev(idx, 0) + dt * v_dfdt(idx, 0);

      // apply Asselin time filter
      if (!leapfrogInit) {
        v_t0(idx, 0) += asselin_ * (v_tp1 - 2.0*v_t0(idx, 0) + v_tm1(idx, 0));
      }

      // move time slices
      v_tm1(idx, 0) = v_t0(idx, 0);
      v_t0(idx, 0) = v_tp1;
    }
  }
  f_t0.set_dirty();
  f_tm1.set_dirty();

  xx.validTime() += tstep_;
}

// -----------------------------------------------------------------------------

void ModelAdvectionBase::finalize(State & xx) const {
  xx_tm1_.clear();
}

// -----------------------------------------------------------------------------

}  // namespace genericMarine
