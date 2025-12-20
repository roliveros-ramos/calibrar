# Adaptative Hierarchical Recombination Evolutionary Strategy (AHR-ES) for derivative-free and black-box optimization

This function performs the optimization of a function using the
Adaptative Hierarchical Recombination Evolutionary Strategy (AHR-ES,
Oliveros & Shin, 2015).

## Usage

``` r
ahres(
  par,
  fn,
  gr = NULL,
  ...,
  lower = -Inf,
  upper = +Inf,
  active = NULL,
  control = list(),
  hessian = FALSE,
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

- active:

  Boolean vector of the same length as par, indicating if the parameter
  is used in the optimization (TRUE) or hold at a fixed value (FALSE).

- control:

  Parameter for the control of the algorithm itself, see details.

- hessian:

  Logical. Should a numerically differentiated Hessian matrix be
  returned? Currently not implemented.

- parallel:

  Logical. Use parallel computation numerical of gradient?

## Value

A list with components:

- par:

  The best set of parameters found.

- value:

  The value of fn corresponding to par.

- counts:

  A two-element integer vector giving the number of calls to fn and gr
  respectively. This excludes those calls needed to compute the Hessian,
  if requested, and any calls to fn to compute a finite-difference
  approximation to the gradient.

- convergence:

  An integer code. 0 indicates successful completion.

- message:

  A character string giving any additional information returned by the
  optimizer, or NULL.

- hessian:

  Only if argument hessian is true. A symmetric matrix giving an
  estimate of the Hessian at the solution found. Note that this is the
  Hessian of the unconstrained problem even if the box constraints are
  active.

## See also

Other optimisers:
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md),
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)

## Author

Ricardo Oliveros-Ramos

## Examples

``` r
if (FALSE) ahres(par=rep(1, 5), fn=sphereN) # \dontrun{}
```
