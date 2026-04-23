# Interpolate time-varying parameters with a spline

Build a smooth time series from values specified at spline knots. This
is useful when a model parameter is allowed to vary through time but is
represented by a smaller set of control values to be estimated.

## Usage

``` r
spline_par(par, n, knots = NULL, periodic = FALSE, period = NULL)

# S3 method for class 'spline.par'
plot(x, ...)
```

## Arguments

- par:

  Numeric vector of parameter values at the knots.

- n:

  Number of output time points.

- knots:

  Optional numeric vector giving the knot positions associated with
  `par`. If `NULL`, equally spaced knots are used.

- periodic:

  Logical. Should the spline be treated as periodic?

- period:

  Optional period of the cycle when `periodic = TRUE` and custom `knots`
  are supplied.

- x:

  A `spline.par` object returned by `spline_par()`.

- ...:

  Further graphical arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

An object of class `spline.par`: a list containing the interpolated
values (`x`), the corresponding time points (`time`), the knot
locations, the knot values, and auxiliary information used by the
plotting method.

## Details

The function rescales knot positions internally to the unit interval
before building the spline. If `periodic = TRUE`, the series is extended
so that the interpolation wraps around smoothly at the boundaries.

The returned object is designed to be plotted with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and to provide
both the discrete interpolated values used in a model and a smoother
curve for visual inspection.

## See also

[`stats::splinefun()`](https://rdrr.io/r/stats/splinefun.html)

## Examples

``` r
sp = spline_par(par = c(0.2, 0.8, 0.4, 0.6), n = 12)
head(sp$x)
#> [1] 0.384375 0.643750 0.778125 0.812500 0.771875 0.681250
plot(sp)

```
