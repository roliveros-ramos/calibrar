library(calibrar)
library(testthat)

tfn = function(x) sum(x^2) + 10
tgr = function(x) 2*x

# optim2 ------------------------------------------------------------------

message("\nTests for optim2() -------- \n")

algorithms = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
               "nlm", "nlminb", "Rcgmin", "Rvmmin", "hjn", 
               "spg", "LBFGSB3", "AHR-ES")

for(alg in algorithms) {
  
  message(alg, " - ", date())
  
  test_that(sprintf("algorithm test: '%s'.", alg), {
    expect_no_error(o0 <- optim2(par=rep(0.5, 2), fn = tfn, method=alg))
    expect_no_error(o1 <- optim2(par=rep(0.5, 2), fn = tfn, gr=tgr, method=alg))
    expect_equal(o0$par, o1$par, tolerance=1e-4)
  })
  
}
