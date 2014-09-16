
# Test functions ----------------------------------------------------------

# method for plotting (2D and 3D)
# method for summary and print
# f() print the minimum

summary.calibrar.function = function(object, ...) {
  return(invisible())
}

SphereN = function(x, sd=0.1, aggregate=TRUE) {
  # f(0,...,0) = 0
  # x_i \in ]-Inf, Inf[
  x = x + rnorm(length(x), sd=sd)
  out = x^2
  if(isTRUE(aggregate)) return(sum(out)) else return(out) 
}
