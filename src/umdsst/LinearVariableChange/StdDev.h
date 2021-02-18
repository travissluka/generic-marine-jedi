/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_LINEARVARIABLECHANGE_STDDEV_H_
#define UMDSST_LINEARVARIABLECHANGE_STDDEV_H_

#include <memory>
#include <string>

#include "oops/util/Printable.h"

// Forward Declaration
namespace eckit {
  class Configuration;
}
namespace umdsst {
  class Geometry;
  class Increment;
  class State;
}

namespace umdsst {

class StdDev: public util::Printable {
 public:
  static const std::string classname() {return "umdsst:StdDev";}

  explicit StdDev(const State &, const State &, const Geometry &,
                  const eckit::Configuration &);

  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  void print(std::ostream &) const override {}
  std::unique_ptr<Increment> stddev_;
};

}  // namespace umdsst

#endif  // UMDSST_LINEARVARIABLECHANGE_STDDEV_H_
