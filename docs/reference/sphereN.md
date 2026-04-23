# Noisy sphere test function

Evaluate a stochastic variant of the sphere test function often used in
optimisation examples. The deterministic sphere function is the sum of
the squared coordinates; here, Gaussian noise is added to each
coordinate before the squared distance is computed.

## Usage

``` r
sphereN(x, sd = 0.1, aggregate = TRUE)
```

## Arguments

- x:

  Numeric vector of coordinates.

- sd:

  Standard deviation of the Gaussian noise added independently to each
  coordinate.

- aggregate:

  Logical. If `TRUE`, return the summed squared distance. If `FALSE`,
  return the contribution of each coordinate separately.

## Value

A numeric scalar if `aggregate = TRUE`, or a numeric vector of squared
contributions if `aggregate = FALSE`.

## Details

This function is mainly useful for examples, testing, and illustrating
how `calibrar` behaves with stochastic objective functions.

## See also

[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md),
[`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md)

## Examples

``` r
set.seed(1)
sphereN(rep(0, 5))
#> [1] 0.03777941
sphereN(rep(0, 5), aggregate = FALSE)
#> [1] 0.0067316837 0.0023758708 0.0054512337 0.0033152416 0.0009326207
```
