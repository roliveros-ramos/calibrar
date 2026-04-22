# Sequential parameter estimation for the calibration of complex models

`calibrate()` minimises an objective function `fn` for the calibration
of complex (possibly expensive and stochastic) models. It supports
sequential calibration in multiple phases (progressively activating
parameters), replicated evaluations for stochastic objectives,
restartable runs, and parallel execution patterns for computationally
intensive models.

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

- phases:

  Optional integer vector of the same length as `par`, indicating the
  phase at which each parameter becomes active. If omitted, all
  parameters are active in a single phase.

- control:

  A list of control options. Common options include `ncores`, `run`,
  `master`, `verbose`, `REPORT`, `restart.file`, `gradient`, and
  `gr.method`. Additional solver-specific options may be passed through
  to the underlying optimiser.

- hessian:

  Logical. Should a numerically differentiated Hessian matrix be
  returned? Currently not implemented.

- replicates:

  Integer or integer vector controlling the number of replicate
  evaluations of `fn` per phase. The default is `1` (deterministic
  objective).

- parallel:

  Logical. Enable parallel computation (e.g., for replicated evaluations
  and/or numerical gradients) using up to `control$ncores` cores.

## Value

An object of class `"calibrar.results"` with components:

- par:

  Best parameter values found (returned with the same structure as input
  `par`).

- value:

  Objective value at `par`.

- counts:

  Number of calls to `fn` and `gr` (where applicable).

- convergence:

  Convergence code returned by the underlying optimiser.

- message:

  Additional information returned by the optimiser, if any.

- method:

  The optimisation method used in the final phase.

- fn:

  The objective function.

- active:

  Logical vector indicating which parameters were active (estimated).

- elapsed:

  Elapsed time for the full calibration run.

- trace:

  Tracing information and per-phase outputs (when available).

## Details

**Sequential phases.** When `phases` is provided, parameters are
activated progressively across phases. In phase \\k\\, parameters with
`phases <= k` are active and estimated, while all others are held fixed.
If `phases` is omitted, all parameters are estimated in a single phase.

**Replicated evaluations.** The argument `replicates` controls the
number of replicate evaluations of `fn` per phase. It can be a single
integer (applied to all phases) or a vector of length equal to the
number of phases. Replication is useful for stochastic models to reduce
Monte Carlo noise and/or to target robust parameter sets.

**Default method selection.** If `method` is not provided, `calibrate()`
defaults to `"Rvmmin"` when `all(replicates == 1)` (deterministic
objective), and to `"AHR-ES"` otherwise (replicated/stochastic
objective).

**Multi-objective outputs.** The objective function `fn` may return a
scalar (single-objective) or a numeric vector of objective components
(multi-objective). Multi-objective optimisation is currently supported
only by methods in `multiMethods` (e.g., `"AHR-ES"`). For
single-objective methods, `fn` must be scalar; alternatively, objective
components can be aggregated into a scalar by defining an aggregated
objective in `calibration_objFn(aggregate = TRUE)`.

**Aggregation.** When `fn` is created via
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
metadata such as the number of components (`nvar`) and component weights
(`weights`) can be stored as attributes. If `aggregate = TRUE`, `fn` is
scalar and can be used with any method. If `aggregate = FALSE` (vector
output), only multi-objective methods can use it directly.

**Parallel execution and run directories.** If `parallel = TRUE`,
`calibrate()` may distribute replicate evaluations and/or
finite-difference gradient computations across multiple cores. The
number of cores is controlled by `control$ncores`. For file-based or
externally executed models, `control$master` and `control$run` can be
used to manage a master/template directory and per-run working
directories (created if needed).

**Restart.** Partial results can be written and used to resume a
calibration via `control$restart.file`. Restart functionality is
currently available only for `"AHR-ES"`, `"Rvmmin"`, and `"hjn"`.

## Choosing an optimisation method

For smooth deterministic objectives, gradient-based methods such as
`"Rvmmin"`, `"L-BFGS-B"`, `"LBFGSB3"`, or `"nlminb"` are often efficient
(with box constraints when needed). For noisy or stochastic objectives
(typically when `replicates > 1`), heuristic/global methods such as
`"AHR-ES"` are generally more appropriate. Not all methods are directly
comparable (local vs.\\ global, deterministic vs.\\ stochastic, scalar
vs.\\ vector-valued objectives); method choice should reflect the
structure of `fn` and the presence of constraints.

## Notes

`"SANN"` is included for compatibility with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), but it is highly
sensitive to tuning and often performs poorly on continuous problems
under default settings. For stochastic objectives, `"AHR-ES"` is the
intended default within `calibrate()`.

## See also

[`calibration_setup`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md),
[`calibration_data`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md),
[`calibration_objFn`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)

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
#> Elapsed time: 0.00s
#> Function value: 0.0258772
#> Parameter values: -1e-04 -1.24e-06 -4.34e-05 3.22e-05 4.01e-05
#> 
#> Status: Rvmminu appears to have converged
#> Optimization using 'Rvmmin' algorithm.
#> Function value: 0.02587717 
#> Status: Rvmminu appears to have converged 
#> Parameters:
#> [1] -9.996345e-05 -1.240568e-06 -4.341452e-05  3.218151e-05  4.013083e-05
#> Computation:
#> function gradient 
#>       54        2 
if (FALSE) { # \dontrun{
calibrate(par=rep(NA, 5), fn=sphereN, replicates=3)
calibrate(par=rep(0.5, 5), fn=sphereN, replicates=3, lower=-5, upper=5)
calibrate(par=rep(0.5, 5), fn=sphereN, replicates=3, lower=-5, upper=5, phases=c(1,1,1,2,3))
calibrate(par=rep(0.5, 5), fn=sphereN, replicates=c(1,1,4), lower=-5, upper=5, phases=c(1,1,1,2,3))
} # }
```
