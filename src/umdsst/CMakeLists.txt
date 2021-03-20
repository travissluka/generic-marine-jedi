#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# the source code in this directory,
# all other files are handled through the subdirectory's respective CMakeLists.txt
list( APPEND umdsst_src_files
  Traits.h
)

# a convenience frunction to make adding source files from subdirectories easier
# (To keep from having to use the full path every time)
function(umdsst_target_sources)
  foreach( src IN LISTS ARGN)
    target_sources(umdsst PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${src})
  endforeach()
endfunction()

find_package( NetCDF REQUIRED COMPONENTS CXX )

# the main library for this interface
ecbuild_add_library( TARGET   umdsst
                     SOURCES  ${umdsst_src_files}
                     LIBS     ${NETCDF_LIBRARIES}
                              eckit eckit_mpi fckit oops ioda ufo saber
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${UMDSST_LINKER_LANGUAGE}
                    )

# add source code in the subdirectories
add_subdirectory(Covariance)
add_subdirectory(Fields)
add_subdirectory(Geometry)
add_subdirectory(GetValues)
add_subdirectory(Increment)
add_subdirectory(LinearVariableChange)
add_subdirectory(ModelAux)
add_subdirectory(State)
add_subdirectory(VariableChange)