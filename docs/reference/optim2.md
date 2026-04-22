# Unified optimisation interface with structured parameters and parallel numerical gradients

`optim2()` provides a unified interface to multiple deterministic and
stochastic optimisation methods, combining the algorithms available
through [`stats::optim()`](https://rdrr.io/r/stats/optim.html) with a
small set of additional solvers accessed via calibrar's internal
dispatcher.

## Usage

``` r
optim2(
  par,
  fn,
  gr = NULL,
  ...,
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "nlm", "nlminb",
    "Rcgmin", "Rvmmin", "hjn", "spg", "LBFGSB3", "AHR-ES"),
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

- method:

  Optimisation method(s) to be used. Can be a single method name or a
  vector of method names (e.g., one per phase). If `NULL`, a default is
  chosen based on `replicates` (see Details).

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

## Details

**Methods.** The current selection includes (i) base
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) methods, (ii)
additional gradient-based solvers from external packages, and (iii)
heuristic/global methods.

**Comparability.** Not all methods are directly comparable: some are
deterministic local optimisers (e.g., quasi-Newton), others are
derivative-free local searches, and others are stochastic/global
heuristics. Method choice should reflect the objective function (smooth
vs. non-smooth, deterministic vs. stochastic) and the presence of
constraints.

**Control arguments.** `optim2()` standardises a set of common control
arguments (e.g., iteration limits, tolerances, tracing, and scaling)
where supported, while still allowing method-specific control parameters
to be passed through to the underlying solver.

**Parallel numerical gradients.** When analytic gradients are not
provided, `optim2()` can compute finite-difference numerical gradients
and distribute function evaluations across multiple cores, which can
substantially reduce wall-clock time for expensive objective functions.

## Choosing an optimisation method

For smooth deterministic objectives, quasi-Newton methods (e.g.,
`"BFGS"` or `"L-BFGS-B"` with box constraints) are often efficient when
gradients are available or can be reliably approximated. For
box-constrained problems, consider `"L-BFGS-B"` or the extended bounded
solvers (e.g., `"Rvmmin"`, `"spg"`). For noisy or stochastic objectives,
heuristic/global methods (e.g., `"AHR-ES"`) may be more appropriate.

## Notes

`"SANN"` is included for compatibility with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), but it is highly
sensitive to tuning and often performs poorly on continuous problems
under default settings. For noisy or rugged objective functions,
`"AHR-ES"` is generally the recommended heuristic alternative within
`optim2()`.

## See also

[`optim`](https://rdrr.io/r/stats/optim.html),
[`nlm`](https://rdrr.io/r/stats/nlm.html),
[`nlminb`](https://rdrr.io/r/stats/nlminb.html)

Other optimisers:
[`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md),
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)

## Author

Ricardo Oliveros-Ramos

## Examples

``` r
optim2(par=rep(NA, 5), fn=sphereN)
#> $par
#> [1] -0.003853539 -0.030986287  0.060781610 -0.040176383  0.060071570
#> 
#> $value
#> [1] 0.003254764
#> 
#> $counts
#> function gradient 
#>      505       NA 
#> 
#> $convergence
#> [1] 1
#> 
#> $message
#> NULL
#> 
```
