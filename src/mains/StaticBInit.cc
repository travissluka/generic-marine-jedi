/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "umdsst/Traits.h"
#include "oops/runs/Run.h"
#include "oops/runs/StaticBInit.h"
#include "saber/oops/instantiateCovarFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<umdsst::Traits>();
  oops::StaticBInit<umdsst::Traits> bmat;
  return run.execute(bmat);
}
