# Sphere function with random noise

This function calculates the Euclidian distance from a point to the
origin after a random displacement of its position.

## Usage

``` r
sphereN(x, sd = 0.1, aggregate = TRUE)
```

## Arguments

- x:

  The coordinates of the point

- sd:

  The standard deviation of the noise to be added to the position of
  `x`, a normal distribution with mean zero is used.

- aggregate:

  If `aggregate` is `TRUE` the distance is returned, otherwise the size
  of the projection of the distance among each axis.

## Value

The distance from the point `x` to the origin after a random
displacement.

## Author

Ricardo Oliveros–Ramos

## Examples

``` r
sphereN(rep(0, 10))
#> [1] 0.03885167
```
