# Objective function between observed and simulated data

Compute an objective value (error, negative log-likelihood, and/or
penalty) given observed data (\`obs\`) and model outputs (\`sim\`). This
function is a thin dispatcher: it resolves \`FUN\` via \[match.fun()\]
and returns \`FUN(obs = obs, sim = sim, ...)\`.

## Usage

``` r
objFn(obs, sim, FUN, ...)

fitness(obs, sim, FUN, ...)
```

## Arguments

- obs:

  Observed data as expected by \`FUN\`. Typically a numeric vector,
  matrix, or array. Missing values (\`NA\`) are permitted; see \*Missing
  values\* below.

- sim:

  Simulated data matching \`obs\`, in the sense expected by \`FUN\`. For
  most pointwise criteria (e.g. \`norm2\`, \`lnorm2\`, \`pois\`),
  \`sim\` should have the same shape as \`obs\`. For composition
  criteria (e.g. \`multinom\`), \`obs\` and \`sim\` are expected to be
  matrices with rows representing samples and columns representing
  classes.

- FUN:

  Objective function to apply. Can be:

  - a function, e.g. \`FUN = norm2\`;

  - or a character string naming a function, e.g. \`FUN = "norm2"\`.

  The function must accept arguments named \`obs\` and \`sim\`
  (additional arguments may be supplied via \`...\`) and must return a
  single numeric scalar.

- ...:

  Additional arguments forwarded to \`FUN\`.

## Value

A numeric scalar: the value of \`FUN(obs = obs, sim = sim, ...)\`. By
convention this is an objective to be \*\*minimised\*\*.

## Details

The returned value is intended to be \*\*minimised\*\* by an optimiser:
lower values indicate a better match between \`sim\` and \`obs\` (or a
lower penalty).

## Conventions and expectations

- \*\*Minimisation\*\*: all provided objectives are formulated so that
  lower values indicate a better fit (or weaker penalty).

- \*\*Shapes\*\*: \`objFn()\` does not reshape or recycle data. It is
  the caller's responsibility to supply \`obs\` and \`sim\` in a
  compatible form for the chosen \`FUN\`.

- \*\*Scalar output\*\*: \`FUN\` should return a length-one numeric
  value. Returning vectors is not supported by \`objFn()\` itself
  (although higher level workflows may aggregate vector-valued
  objectives elsewhere).

## Missing values

Most built-in objectives stop if \`all(is.na(obs))\`. Otherwise, they
ignore missing values using \`na.rm = TRUE\` inside
\`sum()\`/\`mean()\`. This means:

- partial \`NA\`s in \`obs\` are ignored;

- \`NA\`s in \`sim\` will also be dropped from sums when \`na.rm =
  TRUE\` (possibly masking simulation failures if not checked upstream).

## Numerical constraints and stability

Some objectives impose additional constraints:

- \*\*Poisson\*\* (\`pois\`): uses \`log(sim)\`; \`sim\` must be
  strictly positive wherever \`obs \> 0\` (otherwise \`-Inf\` may
  occur). Consider flooring the simulated intensity, e.g. \`sim \<-
  pmax(sim, 1e-12)\`.

- \*\*Log-scale\*\* (\`lnorm2\`, \`lnorm3\`, \`lnorm4\`, \`lnorm4b\`):
  apply \`log(x + tiny)\`; both \`obs + tiny\` and \`sim + tiny\` must
  be positive. If your data may contain zeros, \`tiny\` should be chosen
  accordingly.

- \*\*Compositions\*\* (\`multinom\`): expects matrix inputs and uses
  row sums to convert counts/weights into proportions. See details
  below.

## Supported built-in objective functions

The following functions are available with their current behaviour. All
are formulated as objectives to be minimised.

- `norm2`:

  Sum of squared errors on the original scale: \$\$\sum (obs -
  sim)^2\$\$ Typical use: continuous observations with approximately
  additive errors.

- `lnorm2`:

  Sum of squared errors on the log scale: \$\$\sum (\log(obs + tiny) -
  \log(sim + tiny))^2\$\$ Arguments: \`tiny\` (default \`1e-2\`), added
  before the log to avoid \`log(0)\`. Typical use: positive-valued data
  with multiplicative (lognormal) error structure.

- `lnorm3`:

  Log-scale squared error with an estimated multiplicative scaling
  factor \\q\\. Internally:

  - compute element-wise ratios \`ratio \<- obs/sim\`;

  - set \`NaN\` ratios to \`NA\`;

  - estimate \`q \<- mean(ratio, na.rm = TRUE)\`;

  - return \\\sum (\log(obs+tiny) - \log(sim+tiny) - \log(q))^2\\.

  Typical use: positive data where an overall multiplicative bias (scale
  mismatch) is expected and should not be fully penalised.

- `lnorm4` / `lnorm4b`:

  Extensions of \`lnorm3\` that add a penalty term to discourage extreme
  values of the scaling factor \\q\\ (or extreme per-observation
  ratios). They rely on the helper `rangeq()`:

  - compute \`ratio \<- obs/sim\`, estimate \`q \<- mean(ratio,
    na.rm=TRUE)\`;

  - compute a penalty using parameters \`b\` and \`c\`: \$\$pen = n
    \cdot (\max(\|\log_2(q)\|, b)^c - b^c)\$\$ when \`dump = TRUE\`
    (used by \`lnorm4\`), or \$\$pen = \sum (\max(\|\log_2(ratio)\|,
    b)^c - b^c)\$\$ when \`dump = FALSE\` (used by \`lnorm4b\`), where
    \`n\` is the number of non-missing ratios.

  - add the penalty to the \`lnorm3\`-style objective.

  Arguments: \`tiny\`, \`b\` (default \`1\`), \`c\` (default \`2\`).
  Typical use: log-scale fitting where scale drift must be controlled.

- `pois`:

  Poisson negative log-likelihood (up to constants): \$\$-\sum (obs
  \log(sim) - sim)\$\$ Typical use: counts (or count-like rates) with
  Poisson observation error. Note: \`sim\` must be positive where \`obs
  \> 0\`.

- `multinom`:

  A composition (multinomial-like) objective operating on matrices.
  Inputs are expected as \`obs\` and \`sim\` matrices with:

  - rows = samples (e.g. time steps, hauls, sites),

  - columns = classes (e.g. age/size bins, categories).

  Internal steps (high-level):

  - Let \\A\\ be the number of classes (\`A \<- ncol(sim)\`).

  - Rows of \`sim\` that are all zeros (excluding rows that are all
    \`NA\`) are replaced by \`1\` on that row (interpreted as a uniform
    prior).

  - Row sums are used to compute proportions: \\Psim = sim/sum(sim)\\,
    \\Pobs = obs/sum(obs)\\ (row-wise).

  - Rows with \`sum(sim) == 0\` are set to \`NA\` for numerical
    convenience.

  - Rows with \`sum(obs) == 0\` are set to \`NA\` (interpreted as “no
    proportion data available”).

  - A variance term \`sigma2\` and a small stabiliser \`tiny\` are used
    to define an objective that penalises discrepancies between \`Pobs\`
    and \`Psim\`.

  Arguments: \`size\` (default \`20\`) and \`tiny\` (default \`1e-3\`).
  Interpretation: \`size\` plays the role of an effective sample size
  (larger values typically increase the weight of the composition fit).

- `normp` / `re`:

  Pure penalty on simulated values: \$\$\sum sim^2\$\$ This ignores
  \`obs\` and can be used as a regulariser, or when \`sim\` represents a
  residual vector or deviates already computed upstream. \`re\` is an
  alias of \`normp\`.

- `penalty`:

  Scaled quadratic penalty: \$\$n \cdot mean(sim^2)\$\$ Arguments: \`n\`
  (default \`100\`). This assumes a fixed sample size and can be used to
  put the penalty on a comparable scale across datasets.

## Writing your own objective function

You can supply any custom function via \`FUN\` provided it:

- accepts arguments named \`obs\` and \`sim\` (plus optional \`...\`);

- returns a length-one numeric scalar to be minimised;

- defines its own parameter checks and missing-value policy.

## See also

`match.fun`

## Examples

``` r
## Basic squared-error objective
obs <- c(1, 2, 3, NA, 5)
sim <- c(1.2, 1.9, 2.7, 4.0, 5.1)
objFn(obs, sim, FUN = "norm2")
#> Error in get(as.character(FUN), mode = "function", envir = envir): object 'norm2' of mode 'function' was not found

## Log-scale objective (positive data)
obs <- c(0.1, 1, 10)
sim <- c(0.2, 0.9, 11)
objFn(obs, sim, FUN = lnorm2, tiny = 1e-2)
#> Error: object 'lnorm2' not found

## Poisson objective (counts) with flooring for numerical safety
obs <- c(0, 3, 10, 2)
sim <- c(0, 2.5, 9.8, 1.9)
sim <- pmax(sim, 1e-12)
objFn(obs, sim, FUN = "pois")
#> Error in get(as.character(FUN), mode = "function", envir = envir): object 'pois' of mode 'function' was not found

## Composition objective (matrices: rows = samples, cols = classes)
obs <- rbind(c(10, 5, 0),
            c( 0, 0, 0),  # interpreted as “no composition data”
            c( 2, 1, 7))
sim <- rbind(c( 9, 6, 1),
            c( 0, 0, 0),  # replaced internally by sim+1 on that row
            c( 1, 2, 6))
objFn(obs, sim, FUN = "multinom", size = 20, tiny = 1e-3)
#> Error in get(as.character(FUN), mode = "function", envir = envir): object 'multinom' of mode 'function' was not found

## Custom objective function
my_obj <- function(obs, sim, ...) {
  if (all(is.na(obs))) stop("All observed values are NA.")
  sum(abs(obs - sim), na.rm = TRUE)  # L1 error
}
objFn(obs = c(1, 2, NA), sim = c(1.1, 1.7, 3), FUN = my_obj)
#> [1] 0.4
```
