# spline_par()

#' Interpolate time-varying parameters with a spline
#'
#' Build a smooth time series from values specified at spline knots. This is
#' useful when a model parameter is allowed to vary through time but is
#' represented by a smaller set of control values to be estimated.
#'
#' @param par Numeric vector of parameter values at the knots.
#' @param n Number of output time points.
#' @param knots Optional numeric vector giving the knot positions associated
#'   with `par`. If `NULL`, equally spaced knots are used.
#' @param periodic Logical. Should the spline be treated as periodic?
#' @param period Optional period of the cycle when `periodic = TRUE` and custom
#'   `knots` are supplied.
#'
#' @return
#' An object of class `spline.par`: a list containing the interpolated values
#' (`x`), the corresponding time points (`time`), the knot locations, the knot
#' values, and auxiliary information used by the plotting method.
#'
#' @details
#' The function rescales knot positions internally to the unit interval before
#' building the spline. If `periodic = TRUE`, the series is extended so that the
#' interpolation wraps around smoothly at the boundaries.
#'
#' The returned object is designed to be plotted with `plot()` and to provide
#' both the discrete interpolated values used in a model and a smoother curve
#' for visual inspection.
#'
#' @seealso
#' [stats::splinefun()]
#'
#' @examples
#' sp = spline_par(par = c(0.2, 0.8, 0.4, 0.6), n = 12)
#' head(sp$x)
#' plot(sp)
#'
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

#' @rdname spline_par
#' @param x A `spline.par` object returned by [spline_par()].
#' @param ... Further graphical arguments passed to [graphics::plot()].
#' @export
plot.spline.par = function(x, ...) {
  
  plot(x$plot$time, x$plot$x, type="l", ...)
  points(x$knots, x$par, pch=19)
  if(x$periodic) {
    abline(v=c(0, 1), lty=3, col="blue")
    abline(h=x$plot$zero, lty=2, col="red")
  }
}




# Gaussian dispersion -----------------------------------------------------

# gaussian_kernel()

#' Discretise a bivariate Gaussian kernel on a grid
#'
#' Evaluate a two-dimensional Gaussian density on a regular grid defined by
#' lower and upper bounds in each dimension.
#'
#' @param par A list containing at least `mean` and `sigma`, where `mean` is a
#'   numeric vector of length two and `sigma` is a `2 x 2` covariance matrix.
#' @param lower Numeric vector of length two giving the lower bounds of the
#'   grid.
#' @param upper Numeric vector of length two giving the upper bounds of the
#'   grid.
#' @param n Number of grid cells per dimension. This can be a scalar or a
#'   length-two numeric vector.
#' @param checkSymmetry Logical. Passed to [mvtnorm::dmvnorm()] when checking
#'   the covariance matrix.
#' @param ... Further arguments reserved for future extensions.
#'
#' @return
#' A list with components `x`, `y`, and `z`, suitable for contour or image
#' plotting. `x` and `y` contain the grid coordinates, and `z` contains the
#' corresponding matrix of Gaussian density values.
#'
#' @details
#' This helper is convenient when a smooth bivariate dispersal or movement
#' kernel needs to be represented on a finite grid. The function uses
#' [mvtnorm::dmvnorm()] internally and returns the density values without
#' normalising them to sum to one over the discretised domain.
#'
#' @seealso
#' [mvtnorm::dmvnorm()]
#'
#' @examples
#' ker = gaussian_kernel(
#'   par = list(mean = c(0, 0), sigma = diag(2)),
#'   lower = c(-2, -2),
#'   upper = c(2, 2),
#'   n = c(10, 10)
#' )
#' str(ker)
#'
#' @export
gaussian_kernel = function(par, lower, upper, n=10, checkSymmetry=TRUE, ...) {
  mean  = par$mean
  sigma = par$sigma
  if(length(mean)!=2) stop("The 'mean' vector has to be of length 2.")
  if(nrow(sigma)!=2 | ncol(sigma)!=2) 
    stop("The covariance matrix 'sigma' must have dimension 2.")
  lim = apply(cbind(lower, upper, n), 1, 
              FUN = function(x) pretty(seq(x[1], x[2], length.out=x[3]), x[3]))
  names(lim) =  c("x", "y")
  x = do.call(expand.grid, lim)
  out = dmvnorm(x=x, mean=mean, sigma = sigma, checkSymmetry=checkSymmetry)
  out = array(out, dim=sapply(lim, length))
  out = c(lim, z=list(out))
  return(out)
}

