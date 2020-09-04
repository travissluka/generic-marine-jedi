ecbuild_add_executable( TARGET  oisst_convertstate.x
                        SOURCES ConvertState.cc
                        LIBS    oisst )

ecbuild_add_executable( TARGET  oisst_dirac.x
                        SOURCES Dirac.cc
                        LIBS    oisst )

ecbuild_add_executable( TARGET  oisst_hofx_nomodel.x
                        SOURCES HofXNoModel.cc
                        LIBS    oisst )

ecbuild_add_executable( TARGET  oisst_staticbinit.x
                        SOURCES StaticBInit.cc
                        LIBS    oisst )

ecbuild_add_executable( TARGET  oisst_var.x
                        SOURCES Var.cc
                        LIBS    oisst )
