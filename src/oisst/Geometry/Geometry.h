/*
 * (C) Copyright 2020-2020 UCAR, University of Maryland
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OISST_GEOMETRY_GEOMETRY_H_
#define OISST_GEOMETRY_GEOMETRY_H_

#include <ostream>
#include <string>

#include <memory>

#include "atlas/functionspace.h"
#include "atlas/field.h"

#include "eckit/mpi/Comm.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oisst {
  class GeometryIterator;
}

// ----------------------------------------------------------------------------

namespace oisst {
  // Geometry class
  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
   public:
    static const std::string classname() {return "oisst::Geometry";}

    // constructors and destructor
    explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    // accessors
    const eckit::mpi::Comm & getComm() const {return comm_;}

    // These are needed for the GeometryIterator Interface
    // TODO(template_impl) GeometryIterator begin() const;
    // TODO(template_impl) GeometryIterator end() const;

    atlas::functionspace::StructuredColumns* atlasFunctionSpace() const {
        return atlasFunctionSpace_.get();
    }

   private:
    void print(std::ostream &) const;
    const eckit::mpi::Comm & comm_;

    std::unique_ptr<atlas::functionspace::StructuredColumns>
      atlasFunctionSpace_;
  };
}  // namespace oisst

#endif  // OISST_GEOMETRY_GEOMETRY_H_
