# Numerical gradient approximation

`gradient()` computes a numerical approximation to the gradient of a
scalar objective function at a given parameter vector. It provides
several finite-difference schemes and can distribute the underlying
function evaluations across a registered `foreach` backend when
`parallel = TRUE`.

## Usage

``` r
gradient(fn, x, method, control, parallel, ...)

# Default S3 method
gradient(fn, x, method = NULL, control = list(), parallel = FALSE, ...)
```

## Arguments

- fn:

  Function whose gradient is to be approximated.

- x:

  Numeric vector giving the point at which the gradient is evaluated.

- method:

  Numerical differentiation scheme. If `NULL`, the default method used
  by `gradient.default()` is `"forward"`.

- control:

  A list of control arguments governing step-size construction and, for
  Richardson extrapolation, the extrapolation sequence.

- parallel:

  Logical. Should the function evaluations needed for the numerical
  derivative be carried out in parallel via `foreach`?

- ...:

  Additional arguments passed to `fn`.

## Value

A numeric vector of the same length as `x` containing the estimated
gradient of `fn` at `x`.

## Details

The currently implemented methods are:

- `"forward"` (default): one-sided forward finite differences.

- `"backward"`: one-sided backward finite differences.

- `"central"`: central finite differences.

- `"richardson"`: central differences with Richardson extrapolation.

Step sizes are controlled through `control`. Common entries are:

- `step_method`: how the initial step size is computed. Implemented
  values are `"simple"`, `"optextras"`, and `"start"`.

- `eps`: base perturbation size.

- `zero.tol`: threshold used to decide when a parameter value is treated
  as effectively zero for step-size construction.

- `d`: starting fraction of `x` used by `step_method = "start"`.

- `r` and `v`: Richardson iteration depth and reduction factor for
  `method = "richardson"`.

When `parallel = TRUE`, a suitable `foreach` backend must already be
registered by the user. This function does not itself create or manage a
parallel backend.

## See also

[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md),
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md),
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)

## Examples

``` r
gradient(fn = function(x) sum(x^3), x = c(1, 2), method = "forward")
#> [1]  3 12
gradient(fn = function(x) sum(x^3), x = c(1, 2), method = "central")
#> [1]  3 12
```
