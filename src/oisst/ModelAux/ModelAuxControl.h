/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OISST_MODELAUX_MODELAUXCONTROL_H_
#define OISST_MODELAUX_MODELAUXCONTROL_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oisst {
  class Geometry;
  class ModelAuxIncrement;
}

//-----------------------------------------------------------------------------

namespace oisst {

  // ModelAuxControl class
  class ModelAuxControl : public util::Printable,
                          private util::ObjectCounter<ModelAuxControl> {
   public:
    static const std::string classname() {return "oisst::ModelAuxControl";}

    // constructors / destructor
    ModelAuxControl(const Geometry &, const eckit::Configuration &);
    ModelAuxControl(const Geometry &, const ModelAuxControl &);
    ModelAuxControl(const ModelAuxControl &, const bool);
    ~ModelAuxControl();

    // I/O
    void read(const eckit::Configuration &) {}
    void write(const eckit::Configuration &) const {}

    // math operators
    ModelAuxControl & operator +=(const ModelAuxIncrement &);
    double norm() const {return 0.0;}

   private:
    void print(std::ostream & os) const;
  };
}  // namespace oisst
#endif  // OISST_MODELAUX_MODELAUXCONTROL_H_
