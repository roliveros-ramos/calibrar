# Discretise a bivariate Gaussian kernel on a grid

Evaluate a two-dimensional Gaussian density on a regular grid defined by
lower and upper bounds in each dimension.

## Usage

``` r
gaussian_kernel(par, lower, upper, n = 10, checkSymmetry = TRUE, ...)
```

## Arguments

- par:

  A list containing at least `mean` and `sigma`, where `mean` is a
  numeric vector of length two and `sigma` is a `2 x 2` covariance
  matrix.

- lower:

  Numeric vector of length two giving the lower bounds of the grid.

- upper:

  Numeric vector of length two giving the upper bounds of the grid.

- n:

  Number of grid cells per dimension. This can be a scalar or a
  length-two numeric vector.

- checkSymmetry:

  Logical. Passed to
  [`mvtnorm::dmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)
  when checking the covariance matrix.

- ...:

  Further arguments reserved for future extensions.

## Value

A list with components `x`, `y`, and `z`, suitable for contour or image
plotting. `x` and `y` contain the grid coordinates, and `z` contains the
corresponding matrix of Gaussian density values.

## Details

This helper is convenient when a smooth bivariate dispersal or movement
kernel needs to be represented on a finite grid. The function uses
[`mvtnorm::dmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)
internally and returns the density values without normalising them to
sum to one over the discretised domain.

## See also

[`mvtnorm::dmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)

## Examples

``` r
ker = gaussian_kernel(
  par = list(mean = c(0, 0), sigma = diag(2)),
  lower = c(-2, -2),
  upper = c(2, 2),
  n = c(10, 10)
)
#> Error in do.call(expand.grid, lim): second argument must be a list
str(ker)
#> Error: object 'ker' not found
```
