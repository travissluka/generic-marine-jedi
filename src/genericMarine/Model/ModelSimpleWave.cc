/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Traits.h"
#include "genericMarine/Model/ModelSimpleWave.h"
#include "genericMarine/Geometry/Geometry.h"

#include "oops/interface/ModelBase.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/missingValues.h"

#include "oops/util/abor1_cpp.h"

namespace genericMarine {

// -----------------------------------------------------------------------------

static oops::interface::ModelMaker<Traits, ModelSimpleWave> modelSimpleWave_("SimpleWave");

// -----------------------------------------------------------------------------

ModelSimpleWave::ModelSimpleWave(const Geometry & geom, const ModelSimpleWaveParameters & params)
: geom_(geom), tstep_(params.tstep) {

  // generate the 2D fields that control the advection, this
  // function can be overriden by a subclass. TODO: check to make sure it has the right fields
  params_ = this->setParams();
  if (!params_.has("cx") || !params_.has("cy")) {
    util::abor1_cpp("ModelSimpleWave::setParams must provide \"u\" and \"cy\"",
        __FILE__, __LINE__);
  }
}

// -----------------------------------------------------------------------------

void ModelSimpleWave::print(std::ostream & os) const {
  os << "ModelSimpleWave::print not implemented";
}

// -----------------------------------------------------------------------------

void ModelSimpleWave::initialize(State & xx) const {
  xx_tm1_.clear();
}

// -----------------------------------------------------------------------------

atlas::FieldSet ModelSimpleWave::setParams() const {  
  // This is a simple test situation, of 0 speed at the poles, increasing to ward the equator
  // with a modulation to keep near 0.0 at the coasts
  atlas::FieldSet fset;
  
  // other fields we'll need
  atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());  
  auto v_lonlat = atlas::array::make_view<double, 2>(fspace.lonlat());
  auto v_halo = atlas::array::make_view<int, 1>(fspace.ghost());
  auto v_coastdist = atlas::array::make_view<double, 2>(geom_.extraFields().field("distanceToCoast"));

  // create zero u/v fields
  atlas::Field cx = fspace.createField<double>(atlas::option::name("cx"));
  atlas::Field cy = fspace.createField<double>(atlas::option::name("cy"));
  auto cx_view = atlas::array::make_view<double, 1> (cx);
  auto cy_view = atlas::array::make_view<double, 1> (cy);
  cx_view.assign(0.0);
  cy_view.assign(0.0);
  fset.add(cx);
  fset.add(cy);
    
  // set a horizontally varying u
  const double lat0 = 80.0;
  const double lat1_val = -1.0;
  const double coast_dist = 300e3;
  for(atlas::idx_t idx = 0; idx < fspace.size(); idx++){
    // based on latitude
    double lat = v_lonlat(idx, 1);    
    if (abs(lat) > lat0) {
      cx_view(idx) = 0.0;
    } else {
      cx_view(idx) = (1.0 - abs(lat)/lat0) * lat1_val;
    }

    // set to zero near coast
    cx_view(idx) *= std::min(v_coastdist(idx,0) / coast_dist, 1.0);
  }

  // keep horizontally varying v to 0
  return fset;
}

// -----------------------------------------------------------------------------

void ModelSimpleWave::step(State & xx, const ModelAuxControl &) const {
  atlas::functionspace::StructuredColumns fspace(geom_.functionSpace());  
  double missing; missing = util::missingValue(missing);

  // get various data views we need
  auto dx = atlas::array::make_view<double, 2>(geom_.extraFields().field("dx"));
  auto dy = atlas::array::make_view<double, 2>(geom_.extraFields().field("dy"));
  auto cx = atlas::array::make_view<double, 1>(params_.field("cx"));
  auto cy = atlas::array::make_view<double, 1>(params_.field("cy"));

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
        if (f_tp1(idx,0) == missing) continue;
        
        // dx      
        if (f_t0(idx_xp1, 0) == missing && f_t0(idx_xm1, 0) == missing) {
          // leave as 0, no neighbors
        } else if (f_t0(idx_xm1, 0) == missing) {
          // no left neighbor
        } else if (f_t0(idx_xp1, 0) == missing) {
          // no right neighbor
        } else {
          // 2 neighbors, centered difference          
          dfdx(idx) = (f_t0(idx_xp1, 0) - f_t0(idx_xm1, 0)) / (2.0*dx(idx, 0));
        }


  
  // //         dfdx(idx) = (f_t0(idx_xp1, 0) - f_t0(idx_xm1, 0)) / (2.0*dx(idx, 0));    
  //        }
  //        else {
  //         dfdx(idx) = 3.0;
  //        }
  // //       } else if (f_t0(idx_xp1, 0) == missing) {
  // // //         // dudx = (v_view(idx, 0) - v_view(idx_xm1, 0)) / (2.0*dx_view(idx, 0));
  // //       } else if (f_t0(idx_xm1, 0) == missing) {
  // // //         //dudx = (v_view(idx_xp1, 0) - v_view(idx, 0)) / dx_view(idx, 0);
  // //       }
      }
    }
  }
  
  //--------------------------------------------------
  //TODO I left off here
  //--------------------------------------------------

  // // grid related stuff
  
  // atlas::StructuredGrid grid = fspace.grid();
  // const int size = grid.size();
  // atlas::idx_t ny = grid.ny();
  // double missing;
  // missing = util::missi  {

  // // copy the existing data
  
  // atlas::FieldSet fs_bar = util::copyFieldSet(fs);
  
  // // get views to the fields we'll need  
  // auto c_view = atlas::array::make_view<double,2>(geom_.extraFields().field("phase_speed"));
  // auto lonlat_view = atlas::array::make_view<double, 2>(fspace.lonlat());
  // auto halo_view = atlas::array::make_view<int, 1>(fspace.ghost());
  // auto coastdist_view = atlas::array::make_view<double, 2>(geom_.extraFields().field("distanceToCoast"));

 

  for (int v = 0; v < xx_t0.size(); v++) {
  //   auto v_tm1_view = atlas::array::make_view<double,2>(fs_m1_[v]);
  //   auto v_view = atlas::array::make_view<double,2>(fs[v]);
  //   auto v_bar_view = atlas::array::make_view<double,2>(fs_bar[v]);
    auto v_tp1_view = atlas::array::make_view<double,2>(xx_tp1[v]);
    

  //   // calculate dudx
  //   atlas::Field dudx = fspace.createField<double>();
  //   auto dudx_view = atlas::array::make_view<double, 1>(dudx);
  //   dudx_view.assign(0.0);
    for (atlas::idx_t j = fspace.j_begin(); j < fspace.j_end(); j++) {
      for (atlas::idx_t i = fspace.i_begin(j); i < fspace.i_end(j); i++) {
        atlas::idx_t idx = fspace.index(i,j);        

        //skip land
        if(v_tp1_view(idx,0) == missing) continue;
        v_tp1_view(idx, 0) = dfdx(idx);
  
  //       // horizontal derivative
  //       if (v_view(idx_xp1, 0) == missing && v_view(idx_xm1, 0) == missing) {          
  //       } else if (v_view(idx_xp1, 0) != missing && v_view(idx_xm1, 0) != missing) {          
  //         dudx_view(idx) = (v_view(idx_xp1, 0) - v_view(idx_xm1, 0)) / (2.0*dx_view(idx, 0));    
  //       } else if (v_view(idx_xp1, 0) == missing) {
  //         // dudx = (v_view(idx, 0) - v_view(idx_xm1, 0)) / (2.0*dx_view(idx, 0));
  //       } else if (v_view(idx_xm1, 0) == missing) {
  //         //dudx = (v_view(idx_xp1, 0) - v_view(idx, 0)) / dx_view(idx, 0);
  //       }
      }
    }
  //   fspace.haloExchange(dudx);

  //       // boundary conditions
  //       if (v_view(idx_xp1, 0) == missing && v_view(idx_xm1, 0) != missing) {
  //         dudx_view(idx) = dudx_view(idx_xm1);
  //       } else if (v_view(idx_xp1, 0) != missing && v_view(idx_xm1, 0) == missing) {
  //         dudx_view(idx) = dudx_view(idx_xp1);
  //       }

  //       v_tp1_view(idx, 0) = v_tm1_view(idx, 0) + c * dudx_view(idx) * 2.0*tstep_.toSeconds();
  //       v_bar_view(idx, 0) = v_view(idx, 0);
  //         // + 0.1 * (
  //         //   v_tp1_view(idx, 0) 
  //         // - (2.0*v_view(idx, 0))
  //         // + v_tm1_view(idx, 0));
  //       // v_tp1_view(idx, 0) = dudx_view(idx);
  //     }    
  //   }
     fspace.haloExchange(xx_tp1[v]);
  //   fspace.haloExchange(fs_bar[v]);
  }
  // fs_m1_ = fs_bar;



  xx.validTime() += tstep_;
}

// -----------------------------------------------------------------------------

  void ModelSimpleWave::finalize(State & xx) const {
    xx_tm1_.clear();
  }

  // -----------------------------------------------------------------------------
}  // namespace genericMarine
