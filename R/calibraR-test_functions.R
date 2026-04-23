# Test functions ----------------------------------------------------------

# method for plotting (2D and 3D)
# method for summary and print
# f() print the minimum

#' @export
summary.calibrar.function = function(object, ...) {
  return(invisible())
}

# sphereN()

#' Noisy sphere test function
#'
#' Evaluate a stochastic variant of the sphere test function often used in
#' optimisation examples. The deterministic sphere function is the sum of the
#' squared coordinates; here, Gaussian noise is added to each coordinate before
#' the squared distance is computed.
#'
#' @param x Numeric vector of coordinates.
#' @param sd Standard deviation of the Gaussian noise added independently to
#'   each coordinate.
#' @param aggregate Logical. If `TRUE`, return the summed squared distance. If
#'   `FALSE`, return the contribution of each coordinate separately.
#'
#' @return
#' A numeric scalar if `aggregate = TRUE`, or a numeric vector of squared
#' contributions if `aggregate = FALSE`.
#'
#' @details
#' This function is mainly useful for examples, testing, and illustrating how
#' `calibrar` behaves with stochastic objective functions.
#'
#' @seealso
#' [calibrate()], [optim2()], [ahres()]
#'
#' @examples
#' set.seed(1)
#' sphereN(rep(0, 5))
#' sphereN(rep(0, 5), aggregate = FALSE)
#'
#' @export
sphereN = function(x, sd=0.1, aggregate=TRUE) {
  # f(0,...,0) = 0
  # x_i \in ]-Inf, Inf[
  x = unlist(x)
  x = x + rnorm(length(x), sd=sd)
  out = x^2
  if(isTRUE(aggregate)) return(sum(out)) else return(out) 
}

  
  
  