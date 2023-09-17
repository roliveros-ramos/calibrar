library(calibrar)
library(testthat)

tfn = function(x) sum(x^2) + 10

# optimh ------------------------------------------------------------------

algorithms = c("AHR-ES", "Nelder-Mead", "SANN", "hjn", "LBFGSB3", 
               "cmaes", "genSA", "DE", "soma", "genoud", "PSO", 
               "hybridPSO", "mads", "hjk", "hjkb", "nmk", "nmkb")

for(alg in algorithms) {
  test_that(sprintf("algorithm test: '%s'.", alg), {
    expect_no_error(
      optimh(par=rep(0.5, 5), fn = tfn, lower=rep(-100, 5), upper=rep(100, 5), 
             method=alg))
  })
}

