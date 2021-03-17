ecbuild_add_executable( TARGET  umdsst_convertstate.x
                        SOURCES ConvertState.cc
                        LIBS    umdsst )

ecbuild_add_executable( TARGET  umdsst_dirac.x
                        SOURCES Dirac.cc
                        LIBS    umdsst saber)

ecbuild_add_executable( TARGET  umdsst_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    umdsst )

ecbuild_add_executable( TARGET  umdsst_staticbinit.x
                        SOURCES StaticBInit.cc
                        LIBS    umdsst saber)

ecbuild_add_executable( TARGET  umdsst_var.x
                        SOURCES Var.cc
                        LIBS    umdsst saber)
