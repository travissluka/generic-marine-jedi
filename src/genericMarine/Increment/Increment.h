/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "genericMarine/Fields/Fields.h"

// forward declarations
namespace oops {
  class Variables;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace genericMarine {
  class Geometry;
  class State;
}

// ----------------------------------------------------------------------------

namespace genericMarine {

  // Increment class
  class Increment : private util::ObjectCounter<Increment>,
                    public genericMarine::Fields {
   public:
    static const std::string classname() {return "genericMarine::Increment";}

    // Constructor, destructor
    Increment(const Geometry &, const oops::Variables &,
              const util::DateTime &);
    Increment(const Geometry &, const Increment &);
    Increment(const Increment &, const bool);
    Increment(const Increment &);
    ~Increment();

    // accessors
    std::vector<double> rmsByLevel(const std::string &) const;

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
    void schur_product_with_inv(const Increment &);
    void zero();
    void zero(const util::DateTime &);

    // dirac
    void dirac(const eckit::Configuration &);
  };
}  // namespace genericMarine

