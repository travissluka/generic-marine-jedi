/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OISST_FIELDS_FIELDS_H_
#define OISST_FIELDS_FIELDS_H_

#include <memory>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

// forward declarations
namespace oisst {
  class Geometry;
  class State;
}

namespace atlas {
  class FieldSet;
}

// ----------------------------------------------------------------------------

namespace oisst {

  class Fields : public util::Serializable,
                 public util::Printable {
   public:
    // Constructors/destructors
    Fields(const Geometry & , const oops::Variables & ,
           const util::DateTime &);
    Fields(const Fields &);
    ~Fields() {}

    // math operators
    Fields & operator =(const Fields &);
    Fields & operator+=(const Fields &);
    void accumul(const double &, const Fields &);
    double norm() const;
    void zero();

    // time manipulation
    void updateTime(const util::Duration & dt) { time_ += dt; }
    util::DateTime & validTime() { return time_; }
    const util::DateTime & validTime() const { return time_; }

    // I/O
    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;

    // Serialization (not needed by our project)
    size_t serialSize() const override { return 0; }
    void serialize(std::vector<double> &) const override {}
    void deserialize(const std::vector<double> &, size_t &) override {}

    // other accessors
    std::shared_ptr<atlas::FieldSet> atlasFieldSet() const;
    std::shared_ptr<const Geometry> geometry() const;
    const oops::Variables & variables() const { return vars_; }

   protected:
    std::shared_ptr<atlas::FieldSet> atlasFieldSet_;
    std::shared_ptr<const Geometry> geom_;
    const double missing_;
    util::DateTime time_;
    oops::Variables vars_;

   private:
    void print(std::ostream &) const override;
  };
}  // namespace oisst

#endif  // OISST_FIELDS_FIELDS_H_
