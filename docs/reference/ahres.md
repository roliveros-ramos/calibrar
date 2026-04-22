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

  A numeric vector or list. The length of `par` defines the number of
  parameters to be estimated (i.e., the dimension of the problem).

- fn:

  The objective function to be minimised. It should accept a parameter
  vector (or list, depending on the wrapper) as first argument and
  return either a scalar value (single-objective) or a numeric vector
  (multi-objective).

- gr:

  A function computing the gradient of `fn`. If `NULL`, a numerical
  approximation is used. Alternatively, a character string can specify
  the numerical gradient scheme: `"central"`, `"forward"` (default),
  `"backward"`, or `"richardson"`.

- ...:

  Additional arguments passed to `fn` and `gr`.

- lower:

  Lower bounds for parameters. One value or a vector of the same length
  as `par`. `NA` is treated as `-Inf`. Default is unconstrained.

- upper:

  Upper bounds for parameters. One value or a vector of the same length
  as `par`. `NA` is treated as `Inf`. Default is unconstrained.

- active:

  Boolean vector of the same length as `par`, indicating if the
  parameter is used in the optimisation (`TRUE`) or held at a fixed
  value (`FALSE`).

- control:

  A list of control options. Common options include `ncores`, `run`,
  `master`, `verbose`, `REPORT`, `restart.file`, `gradient`, and
  `gr.method`. Additional solver-specific options may be passed through
  to the underlying optimiser.

- hessian:

  Logical. Should a numerically differentiated Hessian matrix be
  returned? Currently not implemented.

- parallel:

  Logical. Enable parallel computation (e.g., for replicated evaluations
  and/or numerical gradients) using up to `control$ncores` cores.

## Value

A list with components:

- par:

  The best set of parameters found.

- value:

  The value of `fn` corresponding to `par`.

- counts:

  A two-element integer vector giving the number of calls to `fn` and
  `gr` respectively. This excludes those calls needed to compute the
  Hessian, if requested, and any calls to `fn` to compute a
  finite-difference approximation to the gradient.

- convergence:

  An integer code. `0` indicates successful completion.

- message:

  A character string giving any additional information returned by the
  optimizer, or `NULL`.

- hessian:

  Only if argument `hessian` is `TRUE`. A symmetric matrix giving an
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
