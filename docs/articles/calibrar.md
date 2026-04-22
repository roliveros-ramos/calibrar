# Getting started with the \`calibrar\` package

## Introduction

`calibrar` provides general-purpose optimisation tools in R, with an
[`optim()`](https://rdrr.io/r/stats/optim.html)-compatible interface and
additional functionality that is particularly useful when objective
functions are computationally expensive. Although the package was
originally developed for parameter estimation (calibration) of complex
ecological and simulation models, its core optimisation functions can be
used to minimise any user-defined objective function.

The package exposes three complementary entry points:

- [`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md):
  an extension of [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  with a curated set of optimisers and additional features such as
  parameter masking (`active`), optional parallel finite-difference
  gradients, and support for structured (list-based) parameters.
- [`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md):
  a unified interface to a broad collection of heuristic optimisation
  methods, with standardised argument naming.
- [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md):
  a higher-level workflow for sequential (multi-phase) calibration,
  where parameter masking can be progressively relaxed across phases,
  optionally with replicates for stochastic models.

In addition, `calibrar` includes helper tools to build objective
functions from simulated outputs and observational data. These are
optional: advanced users can supply objective functions directly to
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md),
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)
or
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md).

## Basic usage

This vignette introduces three front ends. Choose the one that matches
your task:

- General optimisation with an
  [`optim()`](https://rdrr.io/r/stats/optim.html)-like interface (plus
  masking / parallel finite differences / list parameters): use
  [`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md).
- Heuristic optimisation across many algorithms with standardised
  argument naming: use
  [`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md).
- Sequential (multi-phase) calibration for complex or stochastic models,
  with progressive parameter unmasking: use
  [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md).

### optim2(): general-purpose optimisation (an extension of `stats::optim()`)

[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
uses the same calling convention as
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) for the common
case of a flat numeric parameter vector:

``` r

library(calibrar)

f <- function(x) sum(x^2)

o1 <- stats::optim(par = rep(1, 5), fn = f)
o2 <- optim2(par = rep(1, 5), fn = f)

o1$value; o2$value
#> [1] 1.639553e-07
#> [1] 1.639553e-07
o1$par;   o2$par
#> [1] -1.931714e-04 -3.044063e-04 -1.744066e-04  1.940552e-05  5.641574e-05
#> [1] -1.931714e-04 -3.044063e-04 -1.744066e-04  1.940552e-05  5.641574e-05
```

The results are identical, as here
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
acts just as a wrapper for
[`stats::optim()`](https://rdrr.io/r/stats/optim.html). In addition,
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
can work with structured (list-based) parameters, which can simplify
objective functions for complex models:

``` r

par0 <- list(curve = list(a = 1, b = 0.5), offset = 0)

fn <- function(par) {
(par$curve$a - 2)^2 + (par$curve$b - 1)^2 + (par$offset - 0.1)^2
}

out <- optim2(par = par0, fn = fn)
out$par
#> $curve
#> $curve$a
#> [1] 2.000081
#> 
#> $curve$b
#> [1] 1.00002
#> 
#> 
#> $offset
#> [1] 0.100108
```

The main arguments of
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
are:

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

The first difference is the possible values for the `method` argument.
In addition to the first six methods, also available in
[`optim()`](https://rdrr.io/r/stats/optim.html),
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
gives access to [`stats::nlm()`](https://rdrr.io/r/stats/nlm.html) and
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) but with the
same syntax as [`optim()`](https://rdrr.io/r/stats/optim.html) to make
them easy to use. In addition, three methods from the `optimr` package
(`Rcgmin`, `Rvmmin`, `hjn`), the L-BFGS-B v3 implemented in the
`bfgsb3c` package and the `AHR-ES` (Adaptative Hierarchical
Recombination Evolutionary Strategy) implemented in this package.

We can run the same previous example with two other methods not
available to [`optim()`](https://rdrr.io/r/stats/optim.html):

``` r

optim2(par=rep(1, 5), fn=function(x) sum(x^2), method="nlm")
#> $par
#> [1] -2.500222e-13 -2.500222e-13 -2.500222e-13 -2.500222e-13 -2.500222e-13
#> 
#> $value
#> [1] 3.125556e-25
#> 
#> $counts
#> function gradient 
#>       NA       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] "Relative gradient is close to zero, current iterate is probably solution."
```

``` r

set.seed(880820) # for reproducibility
optim2(par=rep(1, 5), fn=function(x) sum(x^2), method="AHR-ES")
#> $par
#> [1]  3.749723e-10 -2.023475e-10  1.655858e-10 -2.272363e-10  1.488461e-10
#> 
#> $value
#> [1] 2.827589e-19
#> 
#> $counts
#> function gradient 
#>     1024        0 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] "Stopping criteria reached in 128 generations."
#> 
#> $generations
#> [1] 128
```

The second difference is the new argument `active`, which is a vector
indicating if a parameter will be optimized (i.e. active) or fixed to a
constant value during the optimization process. In the next example, we
will fix the third and fourth parameters to its initial values:

``` r

optim2(par=rep(1, 5), fn=function(x) sum(x^2), 
       active=c(TRUE, TRUE, FALSE, FALSE, TRUE))
#> $par
#> [1] -5.962055e-05 -2.965457e-05  1.000000e+00  1.000000e+00 -1.039868e-04
#> 
#> $value
#> [1] 2
#> 
#> $counts
#> function gradient 
#>      126       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
```

As we can see, in the final solution, the `par` value keep the values at
1, the initial value provided. All the numerical gradients computed
internally have also ‘masked’ this parameters and the derivatives are
not computed for them to speed up computation time.

Finally, the third difference is the new argument `parallel`, that
activates the parallel computation of the numerical gradient, when `gr`
is not supplied:

``` r

optim2(par=rep(1, 5), fn=function(x) sum(x^2), parallel=TRUE)
#> $par
#> [1] -1.931714e-04 -3.044063e-04 -1.744066e-04  1.940552e-05  5.641574e-05
#> 
#> $value
#> [1] 1.639553e-07
#> 
#> $counts
#> function gradient 
#>      456       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> NULL
```

This last option will increase performance when the computation time of
`fn` in considerable. We will explain in detail this feature in a
following section.

Additionally, the method for the computation of the numerical gradient
can be chosen within the `control` list:

``` r

optim2(par=rep(0.5, 5), fn=function(x) sum(2*x^(3.1*x)), control=list(gr.method="richardson"))
optim2(par=rep(0.5, 5), fn=function(x) sum(2*x^(3.1*x)), control=list(gr.method="central"))
optim2(par=rep(0.5, 5), fn=function(x) sum(2*x^(3.1*x)), control=list(gr.method="forward"))
```

### optimh()

The function
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)
has a similar functionality as
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
but acts as a wrapper for several heuristic optimization algorithms
implemented in several packages: `dfoptim`, `optimr`, `minqa`, `cmaes`,
`genSA`, `DEoptim`, `soma`, `rgenoud` and `psoptim`. The
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)
function standardizes the inputs and outputs to those of
[`optim()`](https://rdrr.io/r/stats/optim.html), providing a more
convenient user interface. All specific arguments of this methods can be
passed to the original function using the `control` argument.

``` r

optimh(
  par,
  fn,
  gr = NULL,
  ...,
  method = c("AHR-ES", "Nelder-Mead", "SANN", "hjn", "CMA-ES", "genSA", "DE", "soma",
    "genoud", "PSO", "hybridPSO", "mads", "hjk", "hjkb", "nmk", "nmkb"),
  lower = -Inf,
  upper = +Inf,
  active = NULL,
  control = list(),
  hessian = FALSE,
  parallel = FALSE
)
```

``` r

# Covariance Matrix Adaptation Evolutionary Strategy
set.seed(880820) # for reproducibility
optimh(par=rep(1, 5), fn=function(x) sum(x^2), method="CMA-ES",
       control=list(maxit=200))
#> $par
#> [1] -5.817310e-10 -5.739147e-10 -1.575175e-10  3.467173e-10 -1.190319e-10
#> 
#> $value
#> [1] 8.269822e-19
#> 
#> $counts
#> function gradient 
#>     1600       NA 
#> 
#> $convergence
#> [1] 1
#> 
#> $message
#> NULL
```

``` r

# Generalized Simulated Anneling
set.seed(880820) # for reproducibility
optimh(par=rep(1, 5), fn=function(x) sum(x^2), method="genSA", 
       lower=rep(-100, 5), upper=rep(100, 5),
       control=list(maxit=200, temperature=6000))
#> $par
#> [1] 6.988898e-11 6.988898e-11 6.988898e-11 6.988898e-11 6.988898e-11
#> 
#> $value
#> [1] 2.442235e-20
#> 
#> $counts
#> function gradient 
#>     2024        0 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] NA
```

``` r

# Self-Organising Migrating Algorithm
set.seed(880820) # for reproducibility
optimh(par=rep(1, 5), fn=function(x) sum(x^2), method="soma",
       lower=rep(-100, 5), upper=rep(100, 5),
       control=list(maxit=200))
#> $par
#> [1]  1.519459e-16 -4.704646e-16  8.741813e-17  1.169846e-17 -1.104563e-16
#> 
#> $value
#> [1] 2.644038e-31
#> 
#> $counts
#> function gradient 
#>     2000       NA 
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] NA
```

The `maxit` control argument has been standardized to work with all
methods and to represent the maximum number of iterations of the
algorithm. However, the specific number of function evaluations per
iteration may vary between methods. You can refer to the help pages from
each package to have details about every specific method and its control
arguments.

### Running in parallel

Most algorithms implemented in
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
and some in `optimh` can benefit of parallel computation. For the
methods that uses the numerical computation of the gradient, this will
be calculated in parallel. In order to support any type of parallel
implementation, the parallel setup is NOT automatic, and must be done by
the user previous to executed the optimization, as described in the
following example:

``` r

library(parallel)
ncores = detectCores() - 1 # number of cores to be used
cl = makeCluster(ncores)
# this is slower than sequential for very fast models (like this one)
optim2(par=rep(0.5, 5), fn=function(x) sum(x^2), 
               control=list(ncores=ncores), parallel=TRUE)
stopCluster(cl) # close the parallel connections
```

## Sequential parameter estimation

### calibrate()

The
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
function implements an automatic sequential parameter estimation,
meaning parameters can be un-masked (set active) progressively during
sequential phases of the calibration process.

``` r

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

With the basic syntax, the
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
function will work similarly to
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
and
[`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md),
performing a simple optimization:

``` r

calibrate(par=c(1,2,3,NA,NA), fn=function(x) sum(x^2))
#> Using optimization method 'Rvmmin'.
#> Elapsed time: 0.00s
#> Function value: 1.24979e-16
#> Parameter values: -5e-09 -5e-09 -5e-09 -5e-09 -5e-09
#> 
#> Status: Rvmminu appears to have converged
#> Optimization using 'Rvmmin' algorithm.
#> Function value: 1.249787e-16 
#> Status: Rvmminu appears to have converged 
#> Parameters:
#> [1] -4.999645e-09 -4.999291e-09 -4.998935e-09 -5.000000e-09 -5.000000e-09
#> Computation:
#> function gradient 
#>       10        6
```

If `upper` and `lower` bounds are provided, the calibrate function can
take `NA` as starting values for the optimization, to some or all of the
parameters:

``` r

calibrate(par=c(1,2,3,NA,5), fn=function(x) sum(x^2),
          lower=rep(-100, 5), upper=rep(100, 5))
#> Using optimization method 'Rvmmin'.
#> Elapsed time: 0.00s
#> Function value: 1.23526e-16
#> Parameter values: -4.98e-09 -4.97e-09 -4.96e-09 -5e-09 -4.94e-09
#> 
#> Status: Rvmminb appears to have converged
#> Optimization using 'Rvmmin' algorithm.
#> Function value: 1.235263e-16 
#> Status: Rvmminb appears to have converged 
#> Parameters:
#> [1] -4.982589e-09 -4.974351e-09 -4.960025e-09 -4.999959e-09 -4.935026e-09
#> Computation:
#> function gradient 
#>       13        5
```

### Setting up a parameter estimation with multiple phases.

Multiple `phases` can be set up by selecting a different one for each
parameter. The starting value of the optimization for each phase is
updated with the best parameters found in the previous phase:

``` r

calibrate(par=c(1,2,3,NA,5), fn=function(x) sum(x^2),
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,3,2,1))
#> Parameter estimation in three phases.
#> 
#> - Phase 1: 2 out of 5 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 1 finished (0.00s)
#>  Function value: 13
#>  Parameter values: -2.97e-08 -2.61e-08
#> 
#> - Phase 2: 4 out of 5 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 2 finished (0.00s)
#>  Function value: 9
#>  Parameter values: 5.84e-09 -3.93e-08 0 9.41e-09
#> 
#> - Phase 3: 5 out of 5 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 3 finished (0.00s)
#>  Function value: 1.24997e-16
#>  Parameter values: -5e-09 -5e-09 -5e-09 -5e-09 -5e-09
#> 
#> Status: Rvmminb appears to have converged
#> Optimization using 'Rvmmin' algorithm.
#> Function value: 1.24997e-16 
#> Status: Rvmminb appears to have converged 
#> Parameters:
#> [1] -5.000004e-09 -5.000078e-09 -4.999605e-09 -5.000002e-09 -5.000006e-09
#> Computation:
#> function gradient 
#>       11        6
```

When a phase is set to a negative number, the parameter is fixed at its
initial value during all the calibration and it is never optimized:

``` r

calibrate(par=c(1,2,3,NA,5), fn=function(x) sum(x^2),
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,-1,2,1))
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 5 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 1 finished (0.00s)
#>  Function value: 13
#>  Parameter values: -2.97e-08 -2.61e-08
#> 
#> - Phase 2: 4 out of 5 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 2 finished (0.00s)
#>  Function value: 9
#>  Parameter values: 5.84e-09 -3.93e-08 0 9.41e-09
#> 
#> Status: Rvmminb appears to have converged
#> Optimization using 'Rvmmin' algorithm.
#> Function value: 9 
#> Status: Rvmminb appears to have converged 
#> Parameters:
#> [1]  5.835174e-09 -3.934291e-08  3.000000e+00  0.000000e+00  9.410321e-09
#> * Some parameters are not calibrated.
#> Computation:
#> function gradient 
#>        6        4
```

### Dealing with stochastic functions

When dealing with stochastic functions, the argument `replicates` can be
helpful, as it allows to evaluate the objective function several times,
taking the average value of them as the actual function value (as an
approximation of the expected value of the function). When
`replicates=1`, the algorithm used by default is “LBFGSB3”, but when
replicates is greater than 1, “AHR-ES” is used. The next examples use
the function
[`sphereN()`](https://roliveros-ramos.github.io/calibrar/reference/sphereN.md),
which computes the Euclidean distance from a point `x` to the origin of
coordinates after a random displacement of its position (see
[`?sphereN`](https://roliveros-ramos.github.io/calibrar/reference/sphereN.md)
for details). We will set the maximum number of iterations to 1000 to
speed up the execution of the vignette.

``` r

calibrate(par=c(1,2,3,NA,5), fn=sphereN,
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,3,2,1), replicates=3, control=list(maxit=1000))
#> Parameter estimation in three phases.
#> 
#> - Phase 1: 2 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 1 finished (1.27s)
#>  Function value: 12.9126
#>  Parameter values: -0.324 0.12
#> 
#> - Phase 2: 4 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 2 finished (1.25s)
#>  Function value: 10.1169
#>  Parameter values: 0.0695 0.112 0.151 -0.0849
#> 
#> - Phase 3: 5 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 3 finished (1.37s)
#>  Function value: 0.0534613
#>  Parameter values: -0.00843 0.0767 0.0606 0.0265 0.0048
#> 
#> Status: Maximum number of generations or function evaluations reached.
#> Optimization using 'AHR-ES' algorithm.
#> Function value: 0.05346132 
#> Status: Maximum number of generations or function evaluations reached. 
#> Parameters:
#> [1] -0.008425148  0.076711662  0.060580990  0.026532747  0.004796660
#> Computation: (1000 generations)
#> function gradient 
#>     8000        0
```

The number of replicates can be one single value or a vector with the
length equal to the number of phases:

``` r

calibrate(par=c(1,2,3,NA,5), fn=sphereN,
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,3,2,1), replicates=c(1,1,5), control=list(maxit=1000))
#> Parameter estimation in three phases.
#> 
#> - Phase 1: 2 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 1 finished (1.19s)
#>  Function value: 13.3163
#>  Parameter values: -0.346 -0.101
#> 
#> - Phase 2: 4 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 2 finished (1.28s)
#>  Function value: 9.8584
#>  Parameter values: 0.217 -0.335 0.24 0.101
#> 
#> - Phase 3: 5 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 3 finished (1.65s)
#>  Function value: 0.0412229
#>  Parameter values: 0.0282 0.00574 -0.0127 -0.0375 -0.0192
#> 
#> Status: Maximum number of generations or function evaluations reached.
#> Optimization using 'AHR-ES' algorithm.
#> Function value: 0.04122286 
#> Status: Maximum number of generations or function evaluations reached. 
#> Parameters:
#> [1]  0.028236432  0.005735164 -0.012675352 -0.037462605 -0.019210042
#> Computation: (1000 generations)
#> function gradient 
#>     8000        0
```

### Parameters as lists

``` r

calibrate(par=list(par1=c(1,2,3), par2=NA, par3=5), fn=sphereN,
          lower=rep(-100, 5), upper=rep(100, 5),
          phases=c(1,2,-3,2,1), replicates=c(1,5), control=list(maxit=1000))
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 1 finished (1.32s)
#>  Function value: 13.1182
#>  Parameter values: 0.0558 0.181
#> 
#> - Phase 2: 4 out of 5 parameters are currently active.
#>  Using optimization method 'AHR-ES'.
#>  Phase 2 finished (3.22s)
#>  Function value: 8.96767
#>  Parameter values: -0.0167 0.164 -0.14 0.0335
#> 
#> Status: Maximum number of generations or function evaluations reached.
#> Optimization using 'AHR-ES' algorithm.
#> Function value: 8.96767 
#> Status: Maximum number of generations or function evaluations reached. 
#> Parameters:
#>       par11       par12       par13        par2        par3 
#> -0.01674863  0.16445070  3.00000000 -0.13951818  0.03351795 
#> * Some parameters are not calibrated.
#> Computation: (1000 generations)
#> function gradient 
#>     8000        0
```

Note that the function `fn` must be able to take a list as parameter
set, and the user must ensure this works beforehand.

### Running in parallel

Most algorithms implemented in the `calibrate` function can benefit of
parallel computation. For the methods that uses the numerical
computation of the gradient, this will be calculated in parallel as in
[`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md).
In order to support any type of parallel implementation, the parallel
setup is NOT automatic, and must be done by the user previous to
executed the optimization, as described in the following example:

``` r

library(parallel)
ncores = detectCores() - 1 # number of cores to be used
cl = makeCluster(ncores)
# this is slower than sequential for very fast models (like this one)
calib = calibrate(par=rep(0.5, 5), fn=sphereN,
                  replicates=3, 
                  lower=rep(-5, 5), 
                  upper=rep(+5, 5), 
                  phases=c(1,1,1,2,3), 
                  control=list(parallel=TRUE, ncores=ncores))
stopCluster(cl) # close the parallel connections
```

While parallelising an optimization can speed up computations, it is not
always faster than running the optimisation sequentially. Parallel
execution introduces additional overhead, such as data communication,
thread or process management, and synchronization, which can outweigh
the benefits when the objective function is relatively fast to evaluate.
As a result, for computationally inexpensive functions, the time
required to coordinate parallel execution may lead to slower performance
compared to a sequential approach. However, as the evaluation time of
the objective function increases, the efficiency gains from
parallelisation become more significant, making it a valuable strategy
for complex or computationally demanding optimization problems.

Please, refer to `vignette(package="calibrar")` for additional vignettes
or to the [calibrar
website](https://roliveros-ramos.github.io/calibrar/) for more details.
