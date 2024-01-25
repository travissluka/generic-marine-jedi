# (C) COPYRIGHT 2018-2023 UCAR
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

#===================================================================================================
# generic-marine-jedi-bundle
#===================================================================================================

cmake_minimum_required( VERSION 3.12 )
find_package( ecbuild 3.6 REQUIRED)

# project name / version
project( generic-marine-jedi-bundle VERSION 2023.4.0 LANGUAGES C CXX Fortran )

include( ecbuild_bundle )
ecbuild_bundle_initialize()


#===================================================================================================
# required repositories
# JEDI repository versions are from SKYLAB v4 (April 13, 2023)
#===================================================================================================
# Useassume jedi-cmake is in the pre-built environment
include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )

# repos to download
# TODO remove the internal
ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda-internal/oops.git"            TAG 1.9.1 )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda-internal/saber.git"           TAG 1.9.0 )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda-internal/ioda.git"            TAG 2.8.0 )
#ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/jcsda/ioda-converters.git" TAG 0251b27 )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda-internal/ufo.git"             TAG 1.9.0 )
ecbuild_bundle( PROJECT genericMarine   SOURCE ../)

ecbuild_bundle_finalize()
