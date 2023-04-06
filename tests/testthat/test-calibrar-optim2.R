tfn = function(x) sum(x^2) + 10


# optim2 ------------------------------------------------------------------

algorithms = c("Nelder-Mead", "BFGS", "L-BFGS-B", "CG", "SANN",
               "nlm", "nlminb", "Rcgmin", "Rvmmin", "hjn", "spg",
               "LBFGSB3", "AHR-ES")

for(alg in algorithms) {
  test_that(sprintf("optim2 - algorithm test: '%s'.", alg), {
    expect_no_error(
      optim2(par=rep(0.5, 5), fn = tfn, method=alg))
  })
}

