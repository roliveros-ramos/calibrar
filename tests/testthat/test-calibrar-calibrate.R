library(testthat)
library(calibrar)

tfn = function(x) sum(x^2) + 10

foreach::registerDoSEQ()

# calibrate ---------------------------------------------------------------

# algorithms = c('LBFGSB3', 'Rvmmin', 'Rcgmin', 'spg', 'nlminb', 'L-BFGS-B', 
#                'nmkb', 'hjkb', 'mads', 'bobyqa', 'ahr-es', 'cma-es', 'SANN', 
#                'genSA', 'DE', 'soma','genoud', 'PSO', 'PSO2007', 'PSO2011', 
#                'hybridPSO')


algorithms = c("L-BFGS-B", "nlminb", "Rcgmin", "Rvmmin", "hjn", "spg", 
               "LBFGSB3", "cmaes", "genSA", "DE", "soma", "genoud", "PSO", 
               "hybridPSO", "mads", "hjkb", "nmkb", "AHR-ES")

# algorithms = algorithms[1:3]

for(alg in algorithms) {
  
  test_that(sprintf("calibrate - algorithm test: '%s'.", alg), {
    expect_no_error(
      suppressMessages(calibrate(par=rep(0.5, 5), fn = tfn,
                                 lower=rep(-100, 5), upper=rep(100, 5),
                                 method=alg, phases = c(1,1,1,2,1))))
  })
  
  test_that(sprintf("calibrate - algorithm test (parallel): '%s'.", alg), {
    expect_no_error(
      suppressMessages(calibrate(par=rep(0.5, 5), fn = tfn,
                                 lower=rep(-100, 5), upper=rep(100, 5),
                                 method=alg, phases = c(1,1,1,2,1))))
  })
  
}

