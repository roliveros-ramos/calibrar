library(calibrar)
library(testthat)

tfn = function(x) sum(x^2) + 10
tgr = function(x) 2*x

# optim2 ------------------------------------------------------------------

algorithms = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
               "nlm", "nlminb", "Rcgmin", "Rvmmin", "hjn", 
               "spg", "LBFGSB3", "AHR-ES")

for(alg in algorithms) {
  
  test_that(sprintf("algorithm test: '%s'.", alg), {
    expect_no_error(
      optim2(par=rep(0.5, 5), fn = tfn, method=alg))
  })

  test_that(sprintf("algorithm test: '%s'.", alg), {
    expect_no_error(
      optim2(par=rep(0.5, 5), fn = tfn, gr=tgr, method=alg))
  })
  
}
