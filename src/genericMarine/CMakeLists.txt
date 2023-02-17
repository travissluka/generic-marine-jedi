#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# the source code in this directory,
# all other files are handled through the subdirectory's respective CMakeLists.txt
list( APPEND genericmarine_src_files
  Traits.h
)

# a convenience frunction to make adding source files from subdirectories easier
# (To keep from having to use the full path every time)
function(genericmarine_target_sources)
  foreach( src IN LISTS ARGN)
    target_sources(genericmarine PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${src})
  endforeach()
endfunction()

find_package( NetCDF REQUIRED COMPONENTS CXX )

# the main library for this interface
ecbuild_add_library( TARGET   genericmarine
                     SOURCES  ${genericmarine_src_files}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${GENERICMARINE_LINKER_LANGUAGE}
                    )
target_include_directories(genericmarine PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include/genericmarine>)

target_compile_features( genericmarine PUBLIC cxx_std_11 )

target_link_libraries( genericmarine PUBLIC NetCDF::NetCDF_CXX )
target_link_libraries( genericmarine PUBLIC fckit )
target_link_libraries( genericmarine PUBLIC atlas )
target_link_libraries( genericmarine PUBLIC oops )
target_link_libraries( genericmarine PUBLIC saber )
target_link_libraries( genericmarine PUBLIC ioda )
target_link_libraries( genericmarine PUBLIC ufo )

# add source code in the subdirectories
add_subdirectory(Covariance)
add_subdirectory(Fields)
add_subdirectory(Geometry)
add_subdirectory(Increment)
add_subdirectory(LinearVariableChange)
add_subdirectory(ModelAux)
add_subdirectory(State)
add_subdirectory(VariableChange)