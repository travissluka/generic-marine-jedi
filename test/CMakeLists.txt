list( APPEND genericmarine_test_input
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

list( APPEND genericmarine_test_ref
  testref/hofx3d.ref
  testref/dirac.ref
  testref/staticbinit.ref
  testref/var.ref
  )

list( APPEND genericmarine_test_data
  Data/19850101_regridded_sst_1x1.nc
  Data/landmask_1x1.nc
  Data/obs_sst.nc
  Data/rossby_radius.dat
  )


# Create directory for test input/output/reference and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
foreach(FILENAME ${genericmarine_test_input} ${genericmarine_test_ref})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# create static data directory, and symlink all files.
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/Data
          ${CMAKE_CURRENT_BINARY_DIR}/data_static )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated)


# number of PEs to use for MPI enabled tests.
set( GENERICMARINE_TESTS_MAX_MPI 1)


#---------------------------------------------------------------------------------------------------
# The following is a wrapper to simplify the generation of tests of the genericmarine
# executables. There are two types of tests:
#  1) genericmarine executables (EXE must be given)
#  2) interface tests  (SRC must be given)
#
# Arguments:
#  NAME      - the name of the test (after test_genericmarine is prepended), as well
#              as the base name of the input confg file
#  EXE       - name of executable to be used.
#  SRC       - The name of the unit test source file to compile
#  CFG       - The name of the yaml file, if the default testinput/${NAME}.yaml
#              is not to be used
#  MPI       - The number of MPI PEs to use. If not specified, ${GENERICMARINE_TESTS_MAX_MPI}
#              will be used
# TEST_DEPENDS - list of tests this test depends on
#---------------------------------------------------------------------------------------------------
function(genericmarine_add_test)
  # parse the passed arguments
  set(prefix     ARG)
  set(novals     NOCOMPARE NOTRAPFPE)
  set(singlevals NAME EXE SRC CFG MPI CLEAN)
  set(multivals  TEST_DEPENDS)
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

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


  # MPI PEs
  set( MPI ${GENERICMARINE_TESTS_MAX_MPI} )
  if ( ARG_MPI )
    set( MPI ${ARG_MPI})
  endif()

  # # find the MPI command
  # set(MPI_CMD "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_PES}")

  # setup working directory
  set( WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/test_workdir/${ARG_NAME})
  file( MAKE_DIRECTORY ${WORKDIR} )
  foreach(DIRNAME data_generated data_static testinput testoutput testref)
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_BINARY_DIR}/${DIRNAME}
      ${WORKDIR}/${DIRNAME} )
  endforeach()

  # setup output directory
  file( MAKE_DIRECTORY ${WORKDIR}/data_generated/${ARG_NAME})

  # Are we building a unit test / or running a soca executable?
  if ( ARG_SRC )
    # building a unit test
    ecbuild_add_test( TARGET  test_genericmarine_${ARG_NAME}
                      SOURCES executables/${ARG_SRC}
                      ARGS    ${CONFIG_FILE}
                      LIBS    genericmarine
                      MPI     ${MPI}
                      ENVIRONMENT ${TRAPFPE_ENV}
                      WORKING_DIRECTORY ${WORKDIR}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
    set( EXE ./test_genericmarine_${ARG_NAME})

  else()
    # running a genericmarine executable.
    set( EXE ${CMAKE_BINARY_DIR}/bin/${ARG_EXE})
    ecbuild_add_test( TARGET  test_genericmarine_${ARG_NAME}
                      MPI ${MPI}
                      COMMAND ${EXE}
                      WORKING_DIRECTORY ${WORKDIR}
                      ARGS ${CONFIG_FILE}
                      ENVIRONMENT
                              ${TRAPFPE_ENV}
                      DEPENDS ${ARG_EXE}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
  endif()
endfunction()

#================================================================================
# Tests of class interfaces
# TODO: uncomment all of these tests as the appropriate interfaces are implemented
#================================================================================

genericmarine_add_test( NAME geometry
                 SRC TestGeometry.cc )

genericmarine_add_test( NAME state
                 SRC TestState.cc )

genericmarine_add_test( NAME increment
                 SRC TestIncrement.cc )

genericmarine_add_test( NAME modelauxcontrol
                 SRC TestModelAuxControl.cc
                 CFG modelaux.yml )

# genericmarine_add_test( NAME getvalues
#                  SRC TestGetValues.cc )


#    ecbui  ld_add_test(
#      TARGET  test_genericmarine_lineargetvalues
#      SOURCES executables/TestLinearGetValues.cc
#      ARGS    testinput/lineargetvalues.yml
#      MPI     ${MPI_PES}
#      LIBS    genericmarine )

#    ecbuild_add_test(
#     TARGET  test_genericmarine_linearvarchange_stddev
#     SOURCES executables/TestLinearVariableChange.cc
#     ARGS    testinput/linearvarchange_stddev.yml
#     MPI     ${MPI_PES}
#     LIBS    genericmarine )

# ecbuild_add_test(
#   TARGET  test_genericmarine_errorcovariance
#   SOURCES executables/TestErrorCovariance.cc
#   ARGS    testinput/errorcovariance.yml
#   MPI     ${MPI_PES}
#   LIBS    genericmarine
#   TEST_DEPENDS test_genericmarine_staticbinit)

#  ecbuild_add_test(
#    TARGET  test_genericmarine_modelauxcovariance
#    SOURCES executables/TestModelAuxCovariance.cc
#    ARGS    testinput/modelaux.yml
#    MPI     ${MPI_PES}
#    LIBS    genericmarine )

#  ecbuild_add_test(
#    TARGET  test_genericmarine_modelauxincrement
#    SOURCES executables/TestModelAuxIncrement.cc
#    ARGS    testinput/modelaux.yml
#    MPI     ${MPI_PES}
#    LIBS    genericmarine )


#================================================================================
# Test of executables
#================================================================================

genericmarine_add_test( NAME hofx3d
                 EXE  genericmarine_hofx3d.x )

genericmarine_add_test( NAME staticbinit
                 EXE  genericmarine_staticbinit.x )

genericmarine_add_test( NAME dirac
                 EXE  genericmarine_dirac.x
                 TEST_DEPENDS test_genericmarine_staticbinit )

genericmarine_add_test( NAME var
                 EXE  genericmarine_var.x
                 TEST_DEPENDS test_genericmarine_staticbinit )
