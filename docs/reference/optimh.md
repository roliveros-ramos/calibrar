# General-purpose optimisation with heuristic and derivative-free methods

`optimh()` provides a unified interface to the heuristic and
derivative-free optimisation methods available through **calibrar**. It
is intended for objective functions that are expensive, noisy,
non-smooth, stochastic, or otherwise difficult to handle with standard
gradient-based optimisers.

## Usage

``` r
optimh(
  par,
  fn,
  gr = NULL,
  ...,
  method = c("AHR-ES", "Nelder-Mead", "SANN", "hjn", "bobyqa", "CMA-ES", "genSA", "DE",
    "soma", "genoud", "PSO", "hybridPSO", "mads", "hjk", "hjkb", "nmk", "nmkb"),
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

  Function computing the gradient of `fn`. It is included for API
  consistency with
  [`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
  and
  [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
  but is ignored by most methods available through `optimh()`.

- ...:

  Additional arguments passed to `fn` and `gr`.

- method:

  Optimisation method to use. Available choices are `"AHR-ES"`,
  `"Nelder-Mead"`, `"SANN"`, `"hjn"`, `"bobyqa"`, `"CMA-ES"`, `"genSA"`,
  `"DE"`, `"soma"`, `"genoud"`, `"PSO"`, `"hybridPSO"`, `"mads"`,
  `"hjk"`, `"hjkb"`, `"nmk"`, and `"nmkb"`.

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

A list in the same general format as
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md),
typically including:

- par:

  The best parameter values found.

- value:

  The value of `fn` at `par`.

- counts:

  Number of function and gradient evaluations reported by the underlying
  solver, when available.

- convergence:

  Convergence code returned by the underlying solver.

- message:

  Additional information returned by the optimiser, if any.

- hessian:

  A Hessian estimate when supported and requested; otherwise usually
  `NULL`.

## Details

`optimh()` uses the same interface conventions as
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
and
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md):
parameters may be supplied as a vector or list, box constraints can be
enforced through `lower` and `upper`, and subsets of parameters can be
held fixed through `active`.

The available methods include both genuinely heuristic/global algorithms
(for example `"AHR-ES"`, `"CMA-ES"`, `"DE"`, `"genSA"`, `"PSO"`) and a
few derivative-free local searches (for example `"Nelder-Mead"`,
`"bobyqa"`, `"hjn"`, `"nmk"`). The function name should therefore be
read broadly: it is the package entry point for optimisation methods
that are typically useful when classical smooth, gradient-based
optimisation is not the best fit.

Analytic gradients are accepted for interface consistency, but many
methods implemented through `optimh()` ignore `gr`. Likewise, not all
methods return Hessian information even when `hessian = TRUE` is
requested.

## Choosing a method

Heuristic and evolutionary methods are often useful for rugged,
multi-modal, or stochastic objective functions. Derivative-free local
methods may be a better first choice when the objective is deterministic
but gradients are unavailable or unreliable. Method choice should still
reflect the structure of the problem: locality versus global
exploration, box constraints, stochasticity, and computational cost.

## See also

[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md),
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md),
[`gradient()`](https://roliveros-ramos.github.io/calibrar/reference/gradient.md)

Other optimisers:
[`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md),
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)

## Author

Ricardo Oliveros-Ramos

## Examples

``` r
if (FALSE) { # \dontrun{
optimh(
  par = rep(0.5, 5),
  fn = sphereN,
  lower = -5,
  upper = 5,
  method = "AHR-ES",
  control = list(maxit = 200)
)
} # }
```
