ecbuild_add_executable( TARGET  umdsst_convertstate.x
                        SOURCES ConvertState.cc
                        LIBS    umdsst )

ecbuild_add_executable( TARGET  umdsst_dirac.x
                        SOURCES Dirac.cc
                        LIBS    umdsst saber)

ecbuild_add_executable( TARGET  umdsst_hofx_nomodel.x
                        SOURCES HofXNoModel.cc
                        LIBS    umdsst )

ecbuild_add_executable( TARGET  umdsst_staticbinit.x
                        SOURCES StaticBInit.cc
                        LIBS    umdsst saber)

ecbuild_add_executable( TARGET  umdsst_var.x
                        SOURCES Var.cc
                        LIBS    umdsst saber)
