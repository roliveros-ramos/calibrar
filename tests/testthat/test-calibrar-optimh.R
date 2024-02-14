library(calibrar)
library(testthat)

tfn = function(x) sum(x^2) + 10

# optimh ------------------------------------------------------------------

message("\nTests for optimh() -------- \n")

algorithms = c("AHR-ES", "Nelder-Mead", "SANN", "hjn",  
               "CMA-ES", "genSA", "DE", "soma", "genoud", "PSO", 
               "hybridPSO", "mads", "hjk", "hjkb", "nmk", "nmkb", "bobyqa")

for(alg in algorithms) {
  
  message(alg, " - ", date())
  
  test_that(sprintf("algorithm test: '%s'.", alg), {
    expect_no_error(
      optimh(par=rep(0.5, 5), fn = tfn, lower=rep(-100, 5), upper=rep(100, 5), 
             method=alg))
  })
}

message("\nTests for ahres() -------- \n")

test_that("ahres", {
  expect_no_error(
    suppressMessages(ahres(par=rep(0.5, 5), fn = tfn,
                           lower=rep(-100, 5), upper=rep(100, 5)
    )))
})

