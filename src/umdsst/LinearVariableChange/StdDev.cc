/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "umdsst/LinearVariableChange/StdDev.h"
#include "umdsst/State/State.h"
#include "umdsst/Geometry/Geometry.h"
#include "umdsst/Traits.h"
#include "umdsst/Increment/Increment.h"

#include "oops/interface/LinearVariableChange.h"
#include "oops/util/abor1_cpp.h"

namespace umdsst {

// ----------------------------------------------------------------------------

StdDev::StdDev(const State &bkg, const State &traj, const Geometry &geom,
               const eckit::Configuration &conf) {
  stddev_.reset(new Increment(geom, bkg.variables(), bkg.validTime()));
  stddev_->ones();

  if ( conf.has("fixed") ) {
    // a single global fixed value
    double val;
    conf.get("fixed", val);
    *stddev_ *= val;
  } else {
    util::abor1_cpp("StdDev::StdDev() no standard deviation "
                    "method specified", __FILE__, __LINE__);
  }

  // TODO(someone) place the real horizontally varying background error
  // standard deviation in stddev_
}

// ----------------------------------------------------------------------------

void StdDev::multiply(const Increment &dxin, Increment &dxout) const {
  dxout = dxin;
  dxout.schur_product_with(*stddev_);
}

// ----------------------------------------------------------------------------

void StdDev::multiplyAD(const Increment &dxin, Increment &dxout) const {
  multiply(dxin, dxout);
}

// ----------------------------------------------------------------------------

void StdDev::multiplyInverse(const Increment &dxin, Increment &dxout) const {
  dxout = dxin;
  dxout.schur_product_with_inv(*stddev_);
}

// ----------------------------------------------------------------------------

void StdDev::multiplyInverseAD(const Increment &dxin, Increment &dxout) const {
  multiplyInverse(dxin, dxout);
}

// ----------------------------------------------------------------------------

oops::LinearVariableChangeMaker<Traits,
                                oops::LinearVariableChange<Traits, StdDev> >
  makerLinearVariableChangeStdDev_("umdsstStdDev");

// ----------------------------------------------------------------------------

}  // namespace umdsst
