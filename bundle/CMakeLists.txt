# (C) COPYRIGHT 2018-2019 UCAR
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

#
# UMDSST-BUNDLE
#

cmake_minimum_required( VERSION 3.12 )
project( umdsst-bundle VERSION 2021.5.0 LANGUAGES C CXX Fortran )

find_package( ecbuild 3.3.2 QUIET)
include( ecbuild_bundle )
ecbuild_bundle_initialize()


ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda-internal/jedi-cmake.git" UPDATE BRANCH develop )
include( jedicmake/cmake/Functions/git_functions.cmake  )

#===================================================================================================
# ECMWF repos that you probably dont need to build yourself because they should be in the
# jedi-stack and containers
#===================================================================================================
option(BUILD_ECKIT "download and build eckit (not needed if in a jedi container)")
if ( BUILD_ECKIT )
  ecbuild_bundle( PROJECT eckit  GIT "https://github.com/ecmwf/eckit.git"  UPDATE TAG 1.16.0 )
endif ()

option(BUILD_FCKIT "download and build fckit (not needed if in a jedi container)")
if (BUILD_FCKIT)
  ecbuild_bundle( PROJECT fckit  GIT "https://github.com/jcsda/fckit.git"  UPDATE BRANCH release-stable )
  # TODO: replace with https://github.com/ecmwf/fckit.git TAG 0.9.2
endif()

option(BUILD_ATLAS "download and build atlas (not needed if in a jedi container)")
if (BUILD_ATLAS)
  ecbuild_bundle( PROJECT atlas  GIT "https://github.com/jcsda/atlas.git"  UPDATE BRANCH release-stable )
  # TODO: repalce with https://github.com/ecmwf/atlas.git TAG 0.24.1
endif()


#===================================================================================================
# required repositories
#===================================================================================================
ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda/oops.git"            UPDATE BRANCH develop )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda/saber.git"           UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda/ioda.git"            UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/jcsda/ioda-converters.git" UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda/ufo.git"             UPDATE BRANCH develop )
ecbuild_bundle( PROJECT umdsst          SOURCE ../)

ecbuild_bundle_finalize()
