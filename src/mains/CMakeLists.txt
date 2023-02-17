# ecbuild_add_executable( TARGET  genericmarine_convertstate.x
#                         SOURCES ConvertState.cc
#                         LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_staticbinit.x
                        SOURCES StaticBInit.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_dirac.x
                        SOURCES Dirac.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_var.x
                        SOURCES Var.cc
                        LIBS    genericmarine )
