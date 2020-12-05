/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UMDSST_GETVALUES_GEOVALSWRAPPERAD_H_
#define UMDSST_GETVALUES_GEOVALSWRAPPERAD_H_

#include "ufo/GeoVaLs.h"

namespace umdsst {

  extern "C" {
    void geovals_wrapper_ad_fill_f90(const int &, const int &,
                                  const util::DateTime **,
                                  const util::DateTime **,
                                  const atlas::field::FieldImpl *);
  }

  class GeoVaLsWrapperAD {
   public:
    explicit GeoVaLsWrapperAD(const ufo::GeoVaLs & geovals,
                            const ufo::Locations & locs)
     : geovals_(&geovals), locs_(&locs) {}

    void fill(const util::DateTime & t1,
              const util::DateTime & t2,
              const atlas::Field & fld) {
      const util::DateTime * t1p = &t1;
      const util::DateTime * t2p = &t2;
      geovals_wrapper_ad_fill_f90(geovals_->toFortran(), locs_->toFortran(),
                               &t1p, &t2p, fld.get());
    }
   private:
     const ufo::GeoVaLs *geovals_;
     const ufo::Locations *locs_;
  };
}  // namespace umdsst

#endif  // UMDSST_GETVALUES_GEOVALSWRAPPERAD_H_
