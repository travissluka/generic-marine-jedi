# (C) Copyright 2019-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# genericMarine (JEDI based generic DA using file I/O)
################################################################################
cmake_minimum_required( VERSION 3.12 )
project( genericMarine VERSION 2023.2.0 LANGUAGES C CXX Fortran)

find_package( ecbuild 3.6 REQUIRED)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( genericMarine_compiler_flags )


################################################################################
# Dependencies
################################################################################
find_package( NetCDF REQUIRED COMPONENTS Fortran )
find_package( eckit  1.11.6   REQUIRED)
find_package( fckit  0.7.0    REQUIRED)
find_package( atlas  0.20.2   REQUIRED)
find_package( oops   1.5.0    REQUIRED)
find_package( saber  1.5.0    REQUIRED)
find_package( ioda   2.4.0    REQUIRED)
find_package( ufo    1.5.0    REQUIRED)


################################################################################
# Sources
################################################################################
set( GENERICMARINE_LINKER_LANGUAGE CXX )
add_subdirectory( src )
add_subdirectory( test )


################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME genericMarine )

# print the summary of the configuration
ecbuild_print_summary()
