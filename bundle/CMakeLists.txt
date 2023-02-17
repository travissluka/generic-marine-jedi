# (C) COPYRIGHT 2018-2023 UCAR
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

#===================================================================================================
# UMDSST-BUNDLE
#===================================================================================================

cmake_minimum_required( VERSION 3.12 )
find_package( ecbuild 3.6 REQUIRED)

# project name / version
project( genericMarine-bundle VERSION 2023.2.0 LANGUAGES C CXX Fortran )

include( ecbuild_bundle )
ecbuild_bundle_initialize()


#===================================================================================================
# ECMWF repos that you probably dont need to build yourself because they should be in the
# spack-stack and containers
#===================================================================================================
#ecbuild_bundle( PROJECT eckit    GIT "https://github.com/ecmwf/eckit.git" TAG 1.20.2 )
#ecbuild_bundle( PROJECT fckit    GIT "https://github.com/ecmwf/fckit.git" TAG 0.9.5 )
#ecbuild_bundle( PROJECT atlas    GIT "https://github.com/ecmwf/atlas.git" TAG 0.31.1 )

#===================================================================================================
# required repositories
# tags are for JEDI-Skylab v3.0 (Dec 2022)
#===================================================================================================
# Use external jedi-cmake or build in bundle
if(DEFINED ENV{jedi_cmake_ROOT})
  include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )
else()
  ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda/jedi-cmake.git"          UPDATE TAG 1.4.0 )
  include( jedicmake/cmake/Functions/git_functions.cmake )
endif()

ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda/oops.git"            UPDATE TAG 1.5.0 )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda/saber.git"           UPDATE TAG 1.5.0 )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda/ioda.git"            UPDATE TAG 2.4.0 )
#ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/jcsda/ioda-converters.git" UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda/ufo.git"             UPDATE TAG 1.5.0 )
ecbuild_bundle( PROJECT genericMarine   SOURCE ../)

ecbuild_bundle_finalize()
