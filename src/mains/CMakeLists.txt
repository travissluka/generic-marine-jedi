ecbuild_add_executable( TARGET  genericmarine_error_covariance_toolbox.x
                        SOURCES ErrorCovarianceToolbox.cc
                        LIBS    genericmarine)

ecbuild_add_executable( TARGET  genericmarine_forecast.x
                        SOURCES Forecast.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_hofx4d.x
                        SOURCES HofX4D.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_setcorscales.x
                        SOURCES SetCorScales.cc
                        LIBS    genericmarine )

ecbuild_add_executable( TARGET  genericmarine_var.x
                        SOURCES Var.cc
                        LIBS    genericmarine )