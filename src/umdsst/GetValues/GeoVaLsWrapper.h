/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_GETVALUES_GEOVALSWRAPPER_H_
#define UMDSST_GETVALUES_GEOVALSWRAPPER_H_

#include <vector>

#include "ufo/GeoVaLs.h"

#include "oops/util/DateTime.h"
#include "ufo/Locations.h"

namespace umdsst {

extern "C" {
  void geovals_wrapper_fill_f90(const int &,
                                const ufo::Locations &,
                                const util::DateTime **,
                                const util::DateTime **,
                                      atlas::field::FieldImpl **,
                                const int & size);
}

  class GeoVaLsWrapper {
   public:
    explicit GeoVaLsWrapper(const ufo::GeoVaLs & geovals,
                            const ufo::Locations & locs)
     : geovals_(&geovals), locs_(locs) {}

    void fill(const util::DateTime & t1,
              const util::DateTime & t2,
              const std::vector<atlas::Field> & fields) {
      const util::DateTime * t1p = &t1;
      const util::DateTime * t2p = &t2;

      std::vector<atlas::field::FieldImpl*> flds;
      for (int i = 0; i < fields.size(); i++)
        flds.push_back((atlas::field::FieldImpl*)(fields[i].get()));

      geovals_wrapper_fill_f90(geovals_->toFortran(), locs_,
        &t1p, &t2p, flds.data(), flds.size());  // flds.data() <=> ** ;
    }

   private:
     const ufo::GeoVaLs *geovals_;
     const ufo::Locations locs_;
  };
}  // namespace umdsst

#endif  // UMDSST_GETVALUES_GEOVALSWRAPPER_H_
