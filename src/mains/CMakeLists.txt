ecbuild_add_executable( TARGET  genericmarine_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_setcorscales.x
                        SOURCES SetCorScales.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_errorcovariance_training.x
                        SOURCES ErrorCovarianceTraining.cc
                        LIBS    genericmarine)

ecbuild_add_executable( TARGET  genericmarine_dirac.x
                        SOURCES Dirac.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_var.x
                        SOURCES Var.cc
                        LIBS    genericmarine )
