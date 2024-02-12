library(calibrar)
library(testthat)

tfn = function(x) sum(x^3 + 10*x) 

xfn = function(x) x[1]^3 + x[2]^2 + x[1]*x[2] + 3*x[3]
x0 = c(1,1,1)

# gradient ----------------------------------------------------------------

foreach::registerDoSEQ()

methods = c("forward", "backward", "central", "richardson")

for(method in methods) {
  
  test_that(sprintf("gradient - x=1, method test: '%s'.", method), {
    expect_no_error(
      gradient(fn = tfn, x=1, method=method))
  })
  
  test_that(sprintf("gradient (parallel) - x=1,  method test: '%s'.", method), {
    expect_no_error(
      gradient(fn = tfn, x=1, method=method, parallel=TRUE))
  })
  
  test_that(sprintf("gradient - x=[1,1], method test: '%s'.", method), {
    expect_no_error(
      gradient(fn = xfn, x=x0, method=method))
  })
  
  test_that(sprintf("gradient (parallel) - x=[1,1], method test: '%s'.", method), {
    expect_no_error(
      gradient(fn = xfn, x=x0, method=method, parallel=TRUE))
  })
  
  test_that(sprintf("parallel and non-parallel equality, method test: '%s'.", method), {
    expect_equal(gradient(fn = xfn, x=x0, method=method, parallel=TRUE),
                 gradient(fn = xfn, x=x0, method=method))
  })
  
}

test_that(sprintf("gradient - x=1, method test: '%s'.", method), {
  expect_no_error(
    gradient(fn = tfn, x=1, method="central", control=list(step_method="optextras")))
})

