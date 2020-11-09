/*
 * (C) Copyright 2019-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "oisst/Geometry/Geometry.h"
#include "oisst/Increment/Increment.h"
#include "oisst/State/State.h"

#include "eckit/config/Configuration.h"

#include "atlas/field.h"
#include "atlas/array.h"

#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

using atlas::array::make_view;


namespace oisst {

// ----------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom,
                       const oops::Variables & vars,
                       const util::DateTime & vt)
    : Fields(geom, vars, vt) {}

// ----------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom, const Increment & other)
    : Fields(other) {
    // it will normally be used for interpolation and change resolution.
    // For now, just copy without interpolation.
  }

// ----------------------------------------------------------------------------

  Increment::Increment(const Increment & other, const bool copy)
    : Fields(*other.geom_, other.vars_, other.time_) {
    if (copy)
      *this = other;
  }

// ----------------------------------------------------------------------------

  Increment::Increment(const Increment & other)
    : Fields(other) {}

// ----------------------------------------------------------------------------

  Increment::~Increment() { }

// ----------------------------------------------------------------------------

  Increment & Increment::operator =(const Increment &other) {
    // Ligang: Do we assume the geom_ and vars_ are the same, (so they are not
    // empty)? For now yes, as all of the constructors allocate atlasFieldSet.
    // so just assign the values of FieldSet. same for other math operators.
    Fields::operator=(other);
    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator -=(const Increment &other) {
    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
      auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));
      for (int j = 0; j < size; j++)
        fd(j) -= fd_other(j);
    }

    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator +=(const Increment &other) {
    Fields::operator+=(other);
    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator *=(const double &zz) {
    auto fd       = make_view<double, 1>(atlasFieldSet_->field(0));
    const int size = geom_->atlasFunctionSpace()->size();

    for (int j = 0; j < size; j++)
      fd(j) *= zz;

    return *this;
  }

// ----------------------------------------------------------------------------

  void Increment::axpy(const double &zz, const Increment &dx, const bool check)
  {
    ASSERT(!check || time_ == dx.validTime());
    // use accumul, because conceptually it is the same as axpy
    // (axpy is for Increment, accumul is for State)
    accumul(zz, dx);
  }

// ----------------------------------------------------------------------------

  void Increment::diff(const State & x1, const State & x2) {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_x1 = make_view<double, 1>(x1.atlasFieldSet()->field(0));
    auto fd_x2 = make_view<double, 1>(x2.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();

    for (int i = 0; i < size; i++)
      fd(i) = fd_x1(i) - fd_x2(i);
  }

// ----------------------------------------------------------------------------

  double Increment::dot_product_with(const Increment &other) const {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_other = make_view<double, 1>(other.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();
    double dp = 0.0;

    // Ligang: will be updated with missing_value process!
    for (int i = 0; i < size; i++)
      dp += fd(i)*fd_other(i);

    // sum results across PEs
    oops::mpi::world().allReduceInPlace(dp, eckit::mpi::Operation::SUM);

    return dp;
  }

// ----------------------------------------------------------------------------

  void Increment::ones() {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    fd.assign(1.0);
  }

// ----------------------------------------------------------------------------

  void Increment::random() {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    const int size = geom_->atlasFunctionSpace()->size();

    util::NormalDistribution<double> x(size, 0, 1.0, 1);

    for (int i = 0; i < size; i++)
      fd(i) = x[i];
  }

// ----------------------------------------------------------------------------

  void Increment::schur_product_with(const Increment &rhs ) {
    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    auto fd_rhs = make_view<double, 1>(rhs.atlasFieldSet()->field(0));

    const int size = geom_->atlasFunctionSpace()->size();
    for (int i = 0; i < size; i++)
      fd(i) *= fd_rhs(i);
  }

// ----------------------------------------------------------------------------

  void Increment::zero() {
    // Need this wrapper because the overridden zero(time) would otherwise
    // interfere
    Fields::zero();
  }

// ----------------------------------------------------------------------------

  void Increment::zero(const util::DateTime & time) {
    zero();
    time_ = time;
  }

// ----------------------------------------------------------------------------

  void Increment::dirac(const eckit::Configuration & conf) {
    // Ligang: Do NOT know what it is, may ask Jim. Idea of coding from JEDI.
    std::vector<int> ixdir(conf.getIntVector("ixdir"));
    std::vector<int> iydir(conf.getIntVector("iydir"));

    // check
    ASSERT(ixdir.size() > 0 && ixdir.size() == iydir.size());
    const int dir_size = ixdir.size();

    // Ligang: This is where we need the field_data to be 2D;
    // How do we make it 2D? should be related to when create the field.
    const int ny = static_cast<int>(geom_->atlasFunctionSpace()->grid().ny());
    const int nx = static_cast<int>(
      ((atlas::RegularLonLatGrid&)(geom_->atlasFunctionSpace()->grid())).nx() );

    for (int i = 0; i < dir_size; i++)
      ASSERT(ixdir[i] < nx && iydir[i] < ny);

    auto fd = make_view<double, 1>(atlasFieldSet_->field(0));
    for (int i = 0; i < dir_size; i++) {
      int idx = ixdir[i] + iydir[i]*nx;
      fd(idx) = 1.0;
    }
  }

// ----------------------------------------------------------------------------
}  // namespace oisst
