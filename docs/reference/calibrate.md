# Sequential parameter estimation for the calibration of complex models

This function performs the optimization of a function, possibly in
sequential phases of increasing complexity, and it is designed for the
calibration of a model, by minimizing the error function `fn` associated
to it.

## Usage

``` r
calibrate(
  par,
  fn,
  gr,
  ...,
  method,
  lower,
  upper,
  phases,
  control,
  hessian,
  replicates,
  parallel
)

# Default S3 method
calibrate(
  par,
  fn,
  gr = NULL,
  ...,
  method = NULL,
  lower = NULL,
  upper = NULL,
  phases = NULL,
  control = list(),
  hessian = FALSE,
  replicates = 1,
  parallel = FALSE
)
```

## Arguments

- par:

  A numeric vector or list. The length of the par argument defines the
  number of parameters to be estimated (i.e. the dimension of the
  problem).

- fn:

  The function to be minimized.

- gr:

  A function computing the gradient of `fn`. If NULL, a numerical
  approximation of the gradient is used. It can be also a character
  specifying the method for the computation of the numerical gradient:
  'central', 'forward' (the default), 'backward' or 'richardson'.

- ...:

  Additional parameters to be passed to `fn`.

- method:

  The optimization method to be used. The default method is the AHR-ES
  (Adaptative Hierarchical Recombination Evolutionary Strategy,
  Oliveros-Ramos & Shin, 2016). See details for the methods available.

- lower:

  Lower threshold value(s) for parameters. One value or a vector of the
  same length as par. If one value is provided, it is used for all
  parameters. `NA` means `-Inf`. By default `-Inf` is used
  (unconstrained).

- upper:

  Upper threshold value(s) for parameters. One value or a vector of the
  same length as par. If one value is provided, it is used for all
  parameters. `NA` means `Inf`. By default `Inf` is used
  (unconstrained).

- phases:

  An optional vector of the same length as `par`, indicating the phase
  at which each parameter becomes active. If omitted, default value is 1
  for all parameters, performing a single optimization.

- control:

  Parameter for the control of the algorithm itself, see details.

- hessian:

  Logical. Should a numerically differentiated Hessian matrix be
  returned? Currently not implemented.

- replicates:

  The number of replicates for the evaluation of `fn`. The default value
  is 1. A value greater than 1 is only useful for stochastic functions.

- parallel:

  Logical. Use parallel computation numerical of gradient?

## Details

In the control list, `aggFn` is a function to aggregate `fn` to a scalar
value if the returned value is a vector. Some optimization algorithm can
exploite the additional information provided by a vectorial output from
`fn`.

## See also

Other optimisers:
[`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md),
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md),
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)

## Author

Ricardo Oliveros-Ramos

## Examples

``` r
calibrate(par=rep(NA, 5), fn=sphereN)
#> Using optimization method 'Rvmmin'.
#> Elapsed time: 0.05s
#> Function value: 0.0101459
#> Parameter values: 3.93e-05 -2.87e-05 -8.6e-05 -1.04e-05 5.31e-05
#> 
#> Status: Rvmminu appears to have converged
#> Optimization using 'Rvmmin' algorithm.
#> Function value: 0.01014588 
#> Status: Rvmminu appears to have converged 
#> Parameters:
#> [1]  3.927977e-05 -2.873164e-05 -8.598253e-05 -1.039704e-05  5.314853e-05
#> Computation:
#> function gradient 
#>     1057       24 
if (FALSE) { # \dontrun{
calibrate(par=rep(NA, 5), fn=sphereN, replicates=3)
calibrate(par=rep(0.5, 5), fn=sphereN, replicates=3, lower=-5, upper=5)
calibrate(par=rep(0.5, 5), fn=sphereN, replicates=3, lower=-5, upper=5, phases=c(1,1,1,2,3))
calibrate(par=rep(0.5, 5), fn=sphereN, replicates=c(1,1,4), lower=-5, upper=5, phases=c(1,1,1,2,3))
} # }
```
