
# models = c("PoissonMixedModel", "PredatorPrey", "SIR", "IBMLotkaVolterra")
models = c("PoissonMixedModel", "PredatorPrey", "SIR")

tmp = tempdir(check = TRUE)
dir.create(file.path(tmp, "demo"), recursive = TRUE)

message("\nTests for calibrar_demo() -------- \n")

for(model in models) {
  
  message("Testing demo ", model, " - ", date())
  test_that(sprintf("calibrate - algorithm test: '%s'.", model), {
    expect_no_error(suppressMessages(dem <- calibrar_demo(path=file.path(tmp, "demo"), model=model)))
    # expect_no_error(setup <- calibration_setup(file = dem$setup))
    # expect_no_error(observed <- calibration_data(setup=setup, path=dem$path))
  })
  
}

unlink(file.path(tmp, "demo"), recursive = TRUE)    
