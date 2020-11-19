/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_GETVALUES_LOCATIONSWRAPPER_H_
#define UMDSST_GETVALUES_LOCATIONSWRAPPER_H_

#include <memory>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace.h"

#include "ufo/Locations.h"

namespace umdsst {

  extern "C" {
    void locations_wrapper_fill_f90(const int &,
                                    const atlas::field::FieldImpl *);
  }

  class LocationsWrapper {
   public:
    explicit LocationsWrapper(const ufo::Locations & locs) : locs_(locs) {
      atlas::Field field("lonlat",
                         atlas::array::make_datatype<double>(),
                         atlas::array::make_shape(locs.nobs(), 2) );
      locations_wrapper_fill_f90(locs_.toFortran(), field.get());
      functionSpace_.reset(new atlas::functionspace::PointCloud(field));
    }

    const atlas::FunctionSpace & atlasFunctionSpace() const {
        return *functionSpace_;
    }

    const ufo::Locations & locs() const { return locs_; }

   private:
     const ufo::Locations locs_;
     std::unique_ptr<atlas::functionspace::PointCloud> functionSpace_;
  };
}  // namespace umdsst

#endif  // UMDSST_GETVALUES_LOCATIONSWRAPPER_H_
