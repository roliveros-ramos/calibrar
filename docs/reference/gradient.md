# Numerical computation of the gradient, with parallel capabilities

This function calculates the gradient of a function, numerically,
including the possibility of doing it in parallel.

## Usage

``` r
gradient(fn, x, method, control, parallel, ...)
```

## Arguments

- fn:

  The function to calculate the gradient.

- x:

  The value to compute the gradient at.

- method:

  The method used. Currently implemented: central, backward, forward and
  Richardson. See details.

- control:

  A list of control arguments.

- parallel:

  Boolean, should numerical derivatives be calculated in parallel?

- ...:

  Additional arguments to be passed to `fn`.

## Value

The gradient of `fn` at `x`.

## Examples

``` r
gradient(fn=function(x) sum(x^3), x=0)
#> [1] 1e-16
```
