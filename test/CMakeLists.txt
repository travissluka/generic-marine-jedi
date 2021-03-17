list( APPEND umdsst_test_input
  testinput/errorcovariance.yml
  testinput/geometry.yml
  testinput/getvalues.yml
  testinput/hofx3d.yml
  testinput/increment.yml
  testinput/lineargetvalues.yml
  testinput/linearvarchange_stddev.yml
  testinput/modelaux.yml
  testinput/state.yml
  testinput/dirac.yml
  testinput/staticbinit.yml
  testinput/var.yml
  )

list( APPEND umdsst_test_ref
  testref/hofx3d.ref
  testref/dirac.ref
  testref/staticbinit.ref
  testref/var.ref
  )

list( APPEND umdsst_test_data
  Data/19850101_regridded_sst.nc
  Data/obs_sst.nc
  Data/landmask.nc
  )


# create necessary directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/bump)

# link the input files for the tests
foreach(FILENAME
  ${umdsst_test_input}
  ${umdsst_test_ref}
  ${umdsst_test_data}
  )
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()

# number of PEs to use for MPI enabled tests.
set( MPI_PES 2)


#-------------------------------------------------------------------------------
# The following is a wrapper to simplify the generation of tests of the umdsst
# executables.
function(umdsst_exe_test)
  # parse the passed arguments
  set(prefix     ARG)
  set(novals     NOCOMPARE NOTRAPFPE)
  set(singlevals NAME EXE CFG CLEAN)
  set(multivals  TEST_DEPENDS TOL)
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

  # set the tolerances to defaults if none given, or if we want to override
  IF(NOT ARG_TOL)
    SET(ARG_TOL "1.0e-12;0")
  ENDIF()
  LIST(GET ARG_TOL 0 TOL_F)
  LIST(GET ARG_TOL 1 TOL_I)

  # determine if floating point error trapping should be set
  if( ARG_NOTRAPFPE)
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=0")
  else()
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=1")
  endif()

  # determine the default config file name
  if ( ARG_CFG )
    set ( CONFIG_FILE testinput/${ARG_CFG} )
  else()
    set ( CONFIG_FILE testinput/${ARG_NAME}.yml )
  endif()

  # find the MPI command
  set(MPI_CMD "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_PES}")

  # find compare.py provided by oops
  find_file(COMPARE_BIN oops_compare.py HINTS ${oops_BINDIR})

  # add the test
  ecbuild_add_test( TARGET  test_umdsst_${ARG_NAME}
                    TYPE    SCRIPT
                    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/test_wrapper.sh"
                    ARGS    "${CMAKE_BINARY_DIR}/bin/${ARG_EXE}"
                            ${CONFIG_FILE}
                    ENVIRONMENT
                            ${TRAPFPE_ENV}
                            COMPARE_SCRIPT=${COMPARE_BIN}
                            COMPARE_TESTNAME=${ARG_NAME}
                            COMPARE_TOL_F=${TOL_F}
                            COMPARE_TOL_I=${TOL_I}
                            MPI_CMD=${MPI_CMD}
                            SKIP_COMPARE=${ARG_NOCOMPARE}
                            CLEAN_FILES=${ARG_CLEAN}
                    DEPENDS ${ARGE_EXE}
                    TEST_DEPENDS ${ARG_TEST_DEPENDS})
endfunction()

#================================================================================
# Tests of class interfaces
# TODO: uncomment all of these tests as the appropriate interfaces are implemented
#================================================================================

   ecbuild_add_test(
     TARGET  test_umdsst_geometry
     SOURCES executables/TestGeometry.cc
     ARGS    testinput/geometry.yml
     MPI     ${MPI_PES}
     LIBS    umdsst )

   ecbuild_add_test(
     TARGET  test_umdsst_state
     SOURCES executables/TestState.cc
     ARGS    testinput/state.yml
     MPI     ${MPI_PES}
     LIBS    umdsst )

   ecbuild_add_test(
     TARGET  test_umdsst_increment
     SOURCES executables/TestIncrement.cc
     ARGS    testinput/increment.yml
     MPI     ${MPI_PES}
     LIBS    umdsst )

#  ecbuild_add_test(
#    TARGET  test_umdsst_modelauxcontrol
#    SOURCES executables/TestModelAuxControl.cc
#    ARGS    testinput/modelaux.yml
#    MPI     ${MPI_PES}
#    LIBS    umdsst )

   ecbuild_add_test(
     TARGET  test_umdsst_getvalues
     SOURCES executables/TestGetValues.cc
     ARGS    testinput/getvalues.yml
     MPI     ${MPI_PES}
     LIBS    umdsst )

   ecbuild_add_test(
     TARGET  test_umdsst_lineargetvalues
     SOURCES executables/TestLinearGetValues.cc
     ARGS    testinput/lineargetvalues.yml
     MPI     ${MPI_PES}
     LIBS    umdsst )

   ecbuild_add_test(
    TARGET  test_umdsst_linearvarchange_stddev
    SOURCES executables/TestLinearVariableChange.cc
    ARGS    testinput/linearvarchange_stddev.yml
    MPI     ${MPI_PES}
    LIBS    umdsst )

   ecbuild_add_test(
     TARGET  test_umdsst_errorcovariance
     SOURCES executables/TestErrorCovariance.cc
     ARGS    testinput/errorcovariance.yml
     MPI     ${MPI_PES}
     LIBS    umdsst
     TEST_DEPENDS test_umdsst_staticbinit)

#  ecbuild_add_test(
#    TARGET  test_umdsst_modelauxcovariance
#    SOURCES executables/TestModelAuxCovariance.cc
#    ARGS    testinput/modelaux.yml
#    MPI     ${MPI_PES}
#    LIBS    umdsst )

#  ecbuild_add_test(
#    TARGET  test_umdsst_modelauxincrement
#    SOURCES executables/TestModelAuxIncrement.cc
#    ARGS    testinput/modelaux.yml
#    MPI     ${MPI_PES}
#    LIBS    umdsst )


#================================================================================
# Test of executables
#================================================================================

  umdsst_exe_test( NAME hofx3d
                   EXE  umdsst_hofx3d.x )

  umdsst_exe_test( NAME staticbinit
                   EXE  umdsst_staticbinit.x
                   CLEAN Data/bump/* )

  umdsst_exe_test( NAME dirac
                   EXE  umdsst_dirac.x
                   TEST_DEPENDS test_umdsst_staticbinit )

  umdsst_exe_test( NAME var
                   EXE  umdsst_var.x
                   TEST_DEPENDS test_umdsst_staticbinit )
