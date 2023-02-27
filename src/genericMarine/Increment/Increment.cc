/*
 * (C) Copyright 2019-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "genericMarine/Geometry/Geometry.h"
#include "genericMarine/Increment/Increment.h"
#include "genericMarine/State/State.h"

#include "eckit/config/Configuration.h"

#include "atlas/field.h"
#include "atlas/array.h"

#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

using atlas::array::make_view;


namespace genericMarine {

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
    : Fields(other.geom_, other.vars_, other.time_) {
    if (copy) { *this = other; }
  }

// ----------------------------------------------------------------------------

  Increment::Increment(const Increment & other)
    : Fields(other) {}

// ----------------------------------------------------------------------------

  Increment::~Increment() { }

// ----------------------------------------------------------------------------

  std::vector<double> Increment::rmsByLevel(const std::string & varname) const {
    throw eckit::NotImplemented("genericMarine::Increment::rmsByLevel not implemented yet",
                                Here());
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator =(const Increment &other) {
    Fields::operator=(other);
    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator -=(const Increment &other) {
    const int size = geom_.functionSpace().size();

    for (int i = 0; i < vars_.size(); i++) {
      auto fd       = make_view<double, 2>(atlasFieldSet_.field(0));
      auto fd_other = make_view<double, 2>(other.atlasFieldSet_.field(0));
      for (int j = 0; j < size; j++)
        fd(j, 0) -= fd_other(j, 0);
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
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    const int size = geom_.functionSpace().size();

    for (int j = 0; j < size; j++)
      fd(j, 0) *= zz;

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
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    auto fd_x1 = make_view<double, 2>(x1.fieldSet().field(0));
    auto fd_x2 = make_view<double, 2>(x2.fieldSet().field(0));

    const int size = geom_.functionSpace().size();

    for (int i = 0; i < size; i++)
      fd(i, 0) = fd_x1(i, 0) - fd_x2(i, 0);
  }

// ----------------------------------------------------------------------------

  double Increment::dot_product_with(const Increment &other) const {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    auto fd_other = make_view<double, 2>(other.atlasFieldSet_.field(0));
    auto fd_halo = make_view<int, 1>(geom_.functionSpace().ghost());

    const int size = geom_.functionSpace().size();
    double dp = 0.0;

    for (int i = 0; i < size; i++) {
      if (fd_halo(i) != 0) continue;
      dp += fd(i, 0)*fd_other(i, 0);
    }

    // sum results across PEs
    oops::mpi::world().allReduceInPlace(dp, eckit::mpi::Operation::SUM);

    return dp;
  }

// ----------------------------------------------------------------------------

  void Increment::ones() {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    fd.assign(1.0);
  }

// ----------------------------------------------------------------------------

  void Increment::random() {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    auto fd_halo = make_view<int, 1>(geom_.functionSpace().ghost());
    const int size = geom_.functionSpace().size();

    util::NormalDistribution<double> x(size, 0, 1.0, 1);

    // TODO(travis) redo so that random number are computer on PE 0
    // otherwise answers change depending on PE

    int j = 0;
    for (int i = 0; i < size; i++) {
      if (fd_halo(i) != 0) continue;
      fd(i, 0) = x[0];
    }

    atlasFieldSet_.haloExchange();
  }

// ----------------------------------------------------------------------------

  void Increment::schur_product_with(const Increment &rhs ) {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    auto fd_rhs = make_view<double, 2>(rhs.atlasFieldSet_.field(0));

    const int size = geom_.functionSpace().size();
    for (int i = 0; i < size; i++)
      fd(i, 0) *= fd_rhs(i, 0);
  }

// ----------------------------------------------------------------------------

  void Increment::schur_product_with_inv(const Increment &rhs ) {
    auto fd = make_view<double, 2>(atlasFieldSet_.field(0));
    auto fd_rhs = make_view<double, 2>(rhs.atlasFieldSet_.field(0));

    const int size = geom_.functionSpace().size();
    for (int i = 0; i < size; i++)
      fd(i, 0) *= 1.0 / fd_rhs(i, 0);
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
    // get grid related stuff
    const atlas::functionspace::StructuredColumns & fspace =
      static_cast<atlas::functionspace::StructuredColumns>(geom_.functionSpace());
    const int ny = static_cast<int>(fspace.grid().ny());
    const int nx = static_cast<int>(((atlas::RegularLonLatGrid&)(fspace.grid())).nx());
    auto fd_i = make_view<int, 1>(fspace.index_i());
    auto fd_j = make_view<int, 1>(fspace.index_j());
    auto fd_halo = make_view<int, 1>(fspace.ghost());

    // get dirac configuration
    std::vector<eckit::LocalConfiguration> variables;
    conf.get("variables", variables);

    // for each varaible, create diracs
    for ( auto varConf : variables ) {
      // get config for this variable
      std::string varName = varConf.getString("name");
      std::vector<int> ixdir(varConf.getIntVector("ixdir"));
      std::vector<int> iydir(varConf.getIntVector("iydir"));
      const int ndir = ixdir.size();

      // sanity checks
      ASSERT(ixdir.size() > 0 && ixdir.size() == iydir.size());
      for (int i = 0; i < ndir; i++)
        ASSERT(ixdir[i] < nx && iydir[i] < ny);

      // get the relevant field
      auto fd = make_view<double, 2>(atlasFieldSet_.field(varName));

      // for each index
      for (int j = 0; j < ndir; j++) {
        int gidx = -1;
        int ix = ixdir[j]+1;
        int iy = iydir[j]+1;
        for (int i = 0; i < fspace.size(); i++) {
          if (fd_halo(i)) continue;
          if (fd_i(i) == ix && fd_j(i) == iy) {
            fd(i, 0) = 1.0;
            break;
          }
        }
      }
    }
  }

// ----------------------------------------------------------------------------
}  // namespace genericMarine
