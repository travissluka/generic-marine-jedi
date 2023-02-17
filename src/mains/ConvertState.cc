
/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Traits.h"

#include "oops/runs/ConvertState.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::ConvertState<genericMarine::Traits> convertstate;
  return run.execute(convertstate);
}
