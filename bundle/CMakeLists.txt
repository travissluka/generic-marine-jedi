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
project( generic-marine-jedi-bundle VERSION 2023.2.0 LANGUAGES C CXX Fortran )

include( ecbuild_bundle )
ecbuild_bundle_initialize()


#===================================================================================================
# required repositories
# The following tags are are between skylab releases, taken from latest develop dated Feb, 17 2023.
# In the future we will use the skylab tagged releases. (Didn't do it this time because there
# was a very large IODA change I needed to use)
#===================================================================================================
# Useassume jedi-cmake is in the pre-built environment
include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )

# repos to download
ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda/oops.git"            UPDATE TAG 6b49362 )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda/saber.git"           UPDATE TAG 4b1022a )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda/ioda.git"            UPDATE TAG 6f51577 )
#ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/jcsda/ioda-converters.git" UPDATE TAG dfd9c51 )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda/ufo.git"             UPDATE TAG 2e69e83 )
ecbuild_bundle( PROJECT genericMarine   SOURCE ../)

ecbuild_bundle_finalize()
