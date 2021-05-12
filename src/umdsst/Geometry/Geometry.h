/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_GEOMETRY_GEOMETRY_H_
#define UMDSST_GEOMETRY_GEOMETRY_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/field.h"

#include "eckit/mpi/Comm.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace umdsst {
  class GeometryIterator;
}

// ----------------------------------------------------------------------------

namespace umdsst {
  // Geometry class
  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
   public:
    static const std::string classname() {return "umdsst::Geometry";}

    // constructors and destructor
    explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    // load landmask
    void loadLandMask(const eckit::Configuration &);

    // accessors
    const eckit::mpi::Comm & getComm() const {return comm_;}

    // These are needed for the GeometryIterator Interface
    // TODO(template_impl) GeometryIterator begin() const;
    // TODO(template_impl) GeometryIterator end() const;

    atlas::functionspace::StructuredColumns* atlasFunctionSpace() const {
        return atlasFunctionSpace_.get();
    }

    // Ligang: 20210111, adjust for JEDI rep updates.
    atlas::FieldSet* atlasFieldSet() const {
      return atlasFieldSet_.get();
    }

   private:
    atlas::Field interpToGeom(const std::vector<eckit::geometry::Point2> &,
                              const std::vector<double> &) const;
    void readRossbyRadius(const std::string &);
    void print(std::ostream &) const;
    const eckit::mpi::Comm & comm_;

    std::unique_ptr<atlas::functionspace::StructuredColumns>
      atlasFunctionSpace_;
    std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
  };
}  // namespace umdsst

#endif  // UMDSST_GEOMETRY_GEOMETRY_H_
