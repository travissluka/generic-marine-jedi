/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OISST_INCREMENT_INCREMENT_H_
#define OISST_INCREMENT_INCREMENT_H_

#include <string>

#include "oisst/Fields/Fields.h"

#include "oops/base/GeneralizedDepartures.h"

// forward declarations
namespace oops {
  class Variables;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace oisst {
  class Geometry;
  class State;
}

// ----------------------------------------------------------------------------

namespace oisst {

  // Increment class
  class Increment : public oops::GeneralizedDepartures,
                    private util::ObjectCounter<Increment>,
                    public oisst::Fields {
   public:
    static const std::string classname() {return "oisst::Increment";}

    // Constructor, destructor
    Increment(const Geometry &, const oops::Variables &,
              const util::DateTime &);
    Increment(const Geometry &, const Increment &);
    Increment(const Increment &, const bool);
    Increment(const Increment &);
    ~Increment();

    // wrappers of methods that are fully implemented in Fields
    Increment & operator+=(const Increment &);

    // Math operators
    Increment & operator =(const Increment &);
    Increment & operator-=(const Increment &);
    Increment & operator*=(const double &);
    void axpy(const double &, const Increment &, const bool check = true);
    void diff(const State &, const State &);
    double dot_product_with(const Increment &) const;
    void ones();
    void random();
    void schur_product_with(const Increment &);
    void zero();
    void zero(const util::DateTime &);

    // dirac
    void dirac(const eckit::Configuration &);
  };
}  // namespace oisst

#endif  // OISST_INCREMENT_INCREMENT_H_
