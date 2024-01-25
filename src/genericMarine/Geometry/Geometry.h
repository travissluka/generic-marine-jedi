/*
 * (C) Copyright 2020-2022 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/KDTree.h"

#include "eckit/mpi/Comm.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace genericMarine {
  class GeometryIterator;
}

// ----------------------------------------------------------------------------

namespace genericMarine {
  // Geometry class
  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
   public:
    static const std::string classname() {return "genericMarine::Geometry";}

    // constructors and destructor
    explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    // accessors
    const eckit::mpi::Comm & getComm() const {return comm_;}
    void latlon(std::vector<double> &, std::vector<double> &, const bool) const;
    std::vector<size_t> variableSizes(const oops::Variables &) const;
    // These are needed for the GeometryIterator Interface for LETKF
    // TODO(template_impl) GeometryIterator begin() const;
    // TODO(template_impl) GeometryIterator end() const;

    // atlas interfaces
    const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
    atlas::FunctionSpace & functionSpace() {return functionSpace_;}
    const atlas::FieldSet & fields() const {return fields_;}
    atlas::FieldSet & fields() {return fields_;}

    // miscellaneous
    bool levelsAreTopDown() const {return true;}

   private:
    // private helper functions
    atlas::Field interpToGeom(const std::vector<eckit::geometry::Point2> &,
                              const std::vector<double> &) const;
    void loadLandMask(const eckit::Configuration &);
    void readRossbyRadius(const std::string &);
    void calcDistToCoast();
    void print(std::ostream &) const;

    // private variables
    const eckit::mpi::Comm & comm_;
    atlas::functionspace::StructuredColumns functionSpace_;
    atlas::FieldSet fields_;
  };
}  // namespace genericMarine
