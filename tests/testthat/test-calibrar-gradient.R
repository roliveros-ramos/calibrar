tfn = function(x) sum(x^3 + 10*x) 

# gradient ----------------------------------------------------------------

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
      gradient(fn = tfn, x=1, method=method))
  })
  
  test_that(sprintf("gradient (parallel) - x=[1,1], method test: '%s'.", method), {
    expect_no_error(
      gradient(fn = tfn, x=1, method=method, parallel=TRUE))
  })
  
}

