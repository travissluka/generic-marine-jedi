/*
 * (C) Copyright 2020-2023 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

// forward declarations
namespace genericMarine {
  class Geometry;
  class State;
}

// ----------------------------------------------------------------------------

namespace genericMarine {

  class Fields : public util::Serializable,
                 public util::Printable {
   public:
    // Constructors/destructors
    Fields(const Geometry & , const oops::Variables & ,
           const util::DateTime &);
    Fields(const Fields &);
    ~Fields() {}
    void updateFields(const oops::Variables &);

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
    // std::shared_ptr<atlas::FieldSet> atlasFieldSet() const;
    // std::shared_ptr<const Geometry> geometry() const;
    const oops::Variables & variables() const { return vars_; }
    const Geometry & geometry() const { return geom_; }
    const util::DateTime & time() const { return time_; }

    // atlas interfaces
    void toFieldSet(atlas::FieldSet &) const;
    void toFieldSetAD(const atlas::FieldSet &);
    void fromFieldSet(const atlas::FieldSet &);

    atlas::FieldSet & fieldSet() { return atlasFieldSet_;}
    const atlas::FieldSet & fieldSet() const { return atlasFieldSet_;}

   protected:
    atlas::FieldSet atlasFieldSet_;
    const Geometry & geom_;
    const double missing_;
    util::DateTime time_;
    oops::Variables vars_;

   private:
    void print(std::ostream &) const override;
  };
}  // namespace genericMarine
