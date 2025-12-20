# Predict time-varying parameters using splines.

Predict time-varying parameters using splines.

## Usage

``` r
spline_par(par, n, knots = NULL, periodic = FALSE, period = NULL)
```

## Arguments

- par:

  Values at knots

- n:

  Number of points. Time (independent variable) is assumed to be between
  0 and n with length(par) equidistant points (including 0 and n).

- knots:

  Position of knots. Default, is length(x) equidistant points between 0
  and 1. Always are re-scaled to 0 to 1.

- periodic:

  boolean, is the spline periodic?

- period:

  If periodic is TRUE, it specify the time period.

## Value

A list with the interpolates values as 'x' and 'time'.
