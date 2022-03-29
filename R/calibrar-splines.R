#' Predict time-varying parameters using splines.
#'
#' @param x Values at knots
#' @param n Number of knots
#' @param periodic boolean, is the spline periodic?
#'
#' @return
#' @export
spline_par = function(par, n, periodic=FALSE) {

  if(isTRUE(periodic)) par = c(par, par[1])
  
  it = seq(from=0, to=n, length.out=length(par))
  pt = seq(from=0.5, to=n)
  
  spfun = splinefun(x=it, y=par)
  out = spfun(pt)

  out = list(time=pt, x=out, otime=it, par=par)
  class(out) = "spline.par"
  
  return(out)
  
}


#' @export
plot.spline.par = function(x, ...) {
  
  plot(x$time, x$x, type="l", ...)
  points(x$otime, x$par, pch=19)
  
}

