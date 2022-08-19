#' Predict time-varying parameters using splines.
#'
#' @param x Values at knots
#' @param n Number of points. Time (independent variable) is assumed to be between 0 and n with length(par) equidistant points (including 0 and n).
#' @param knots Position of knots. Default, is length(x) equidistant points between 0 and 1. Always are re-scaled to 0 to 1.
#' @param periodic boolean, is the spline periodic?
#'
#' @return A list with the interpolates values as 'x' and 'time'.
#' @export
spline_par = function(par, n, knots=NULL, periodic=FALSE) {

  if(isTRUE(periodic)) par = c(par, par[1])
  if(!is.null(knots)) {
    if(any(is.na(knots))) stop("'knots' may not include NAs.")
    if(length(knots)!=length(par)) stop("'knots' and 'par' must have the same length.")
    if(any(diff(knots)<=0)) stop("'knots' must be strictly increasing.")
    knots = (knots - min(knots))/(max(knots)-min(knots))
  } else {
    knots = seq(from=0, to=n, length.out=length(par))/n
  }
  
  pt = seq(from=0.5, to=n, by=1)/n
  
  spfun = splinefun(x=knots, y=par)
  out   = spfun(pt)

  out = list(time=pt, x=out, knots=knots, par=par)
  class(out) = "spline.par"
  
  return(out)
  
}

#' @export
plot.spline.par = function(x, ...) {
  
  plot(x$time, x$x, type="l", ...)
  points(x$knots, x$par, pch=19)
  
}

