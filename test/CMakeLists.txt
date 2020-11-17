list( APPEND umdsst_test_input
  testinput/geometry.yml
  testinput/errorcovariance.yml
  testinput/getvalues.yml
  testinput/lineargetvalues.yml
  testinput/linearmodel.yml
  testinput/increment.yml
  testinput/model.yml
  testinput/modelaux.yml
  testinput/state.yml
  )


# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${umdsst_test_input})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()


list( APPEND umdsst_test_data
  Data/19850101_regridded_sst.nc
  Data/obs_sst.nc
  )


# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
foreach(FILENAME ${umdsst_test_data})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()

# number of PEs to use for MPI enabled tests.
set( MPI_PES 2)

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

#  ecbuild_add_test(
#    TARGET  test_umdsst_getvalues
#    SOURCES executables/TestGetValues.cc
#    ARGS    testinput/getvalues.yml
#    MPI     ${MPI_PES}
#    LIBS    umdsst )

#  ecbuild_add_test(
#    TARGET  test_umdsst_lineargetvalues
#    SOURCES executables/TestLinearGetValues.cc
#    ARGS    testinput/lineargetvalues.yml
#    MPI     ${MPI_PES}
#    LIBS    umdsst )

#  ecbuild_add_test(
#    TARGET  test_umdsst_errorcovariance
#    SOURCES executables/TestErrorCovariance.cc
#    ARGS    testinput/errorcovariance.yml
#    MPI     ${MPI_PES}
#    LIBS    umdsst )

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
