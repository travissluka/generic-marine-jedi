# (C) COPYRIGHT 2018-2025 UCAR
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

#===================================================================================================
# generic-marine-jedi-bundle
#===================================================================================================

cmake_minimum_required( VERSION 3.12 )
find_package( ecbuild 3.6 REQUIRED)

# project name / version
project( generic-marine-jedi-bundle VERSION 2025.1.0 LANGUAGES C CXX Fortran )

include( ecbuild_bundle )
ecbuild_bundle_initialize()


#===================================================================================================
# required repositories
# JEDI repository versions from develop branches on 10/21/2024
#===================================================================================================
# Useassume jedi-cmake is in the pre-built environment
include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )

# repos to download
# TODO remove the internal
ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda/oops.git"            TAG 81410ff )
ecbuild_bundle( PROJECT vader           GIT "https://github.com/jcsda/vader.git"           TAG aedb250 )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda/saber.git"           TAG d60ee2c )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda/ioda.git"            TAG 31dc28f )
#ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/jcsda/ioda-converters.git" TAG 0251b27 ) # TODO update this
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda/ufo.git"             TAG a3b705a )
ecbuild_bundle( PROJECT genericMarine   SOURCE ../)

ecbuild_bundle_finalize()
