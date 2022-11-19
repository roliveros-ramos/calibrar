#' Predict time-varying parameters using splines.
#'
#' @param x Values at knots
#' @param n Number of points. Time (independent variable) is assumed to be between 0 and n with length(par) equidistant points (including 0 and n).
#' @param knots Position of knots. Default, is length(x) equidistant points between 0 and 1. Always are re-scaled to 0 to 1.
#' @param periodic boolean, is the spline periodic?
#'
#' @return A list with the interpolates values as 'x' and 'time'.
#' @export
spline_par = function(par, n, knots=NULL, periodic=FALSE, period=NULL) {

  if(!is.null(knots)) {
    if(any(is.na(knots))) stop("'knots' may not include NAs.")
    if(length(knots)!=length(par)) stop("'knots' and 'par' must have the same length.")
    if(any(diff(knots)<=0)) stop("'knots' must be strictly increasing.")
    if(isTRUE(periodic)) {
      if(is.null(period)) stop("A value for 'period' must be indicated if knots are provided")
      knots = (knots/period) %% period
      if(any(duplicated(knots))) stop("Knots must not be duplicated, check values after periodicity taken into account.")
      ind = order(knots)
      knots = knots[ind]
      par   = par[ind]
      if(all(c(0, 1) %in% knots)) stop("Redundant values for extremes in periodic spline has been supplied.")
      oknots = knots
      opar = par
      knots = c(knots-1, knots, knots+1)
      par = c(par, par, par)
    } else {
      knots = (knots - min(knots))/(max(knots)-min(knots))
    }
  } else {
    if(isTRUE(periodic)) {
      par = c(par, par[1])
    }
    knots = seq(from=0, to=n, length.out=length(par))/n
  }
  
  pt = seq(from=0.5, to=n, by=1)/n
  xpt = seq(from=-0.5, to=n+0.5, by=0.5)/n
  
  spfun = splinefun(x=knots, y=par)
  out   = spfun(pt)
  xout = spfun(xpt)
  pval = spfun(0)

  if(isTRUE(periodic)) {
    par = opar
    knots = oknots
  }
  
  out = list(time=pt, x=out, knots=knots, par=par,
             plot=list(time=xpt, x=xout, zero=pval), periodic=periodic)
  class(out) = "spline.par"
  
  return(out)
  
}

#' @export
plot.spline.par = function(x, ...) {
  
  plot(x$plot$time, x$plot$x, type="l", ...)
  points(x$knots, x$par, pch=19)
  if(x$periodic) {
    abline(v=c(0, 1), lty=3, col="blue")
    abline(h=x$plot$zero, lty=2, col="red")
  }
}

