list( APPEND oisst_test_input
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
foreach(FILENAME ${oisst_test_input})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()


#================================================================================
# Tests of class interfaces
# TODO: uncomment all of these tests as the appropriate interfaces are implemented
#================================================================================

# ecbuild_add_test(
#   TARGET  test_oisst_geometry
#   SOURCES executables/TestGeometry.cc
#   ARGS    testinput/geometry.yml
#   LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_state
#    SOURCES executables/TestState.cc
#    ARGS    testinput/state.yml
#    LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_increment
#    SOURCES executables/TestIncrement.cc
#    ARGS    testinput/increment.yml
#    LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_modelauxcontrol
#    SOURCES executables/TestModelAuxControl.cc
#    ARGS    testinput/modelaux.yml
#    LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_getvalues
#    SOURCES executables/TestGetValues.cc
#    ARGS    testinput/getvalues.yml
#    LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_lineargetvalues
#    SOURCES executables/TestLinearGetValues.cc
#    ARGS    testinput/lineargetvalues.yml
#    LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_errorcovariance
#    SOURCES executables/TestErrorCovariance.cc
#    ARGS    testinput/errorcovariance.yml
#    LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_modelauxcovariance
#    SOURCES executables/TestModelAuxCovariance.cc
#    ARGS    testinput/modelaux.yml
#    LIBS    oisst )

#  ecbuild_add_test(
#    TARGET  test_oisst_modelauxincrement
#    SOURCES executables/TestModelAuxIncrement.cc
#    ARGS    testinput/modelaux.yml
#    LIBS    oisst )
