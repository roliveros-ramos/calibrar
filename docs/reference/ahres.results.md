# Evolutionary-strategy result objects

Objects of class `ahres.results` are returned by
[`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md).
They store the fitted parameters and the main diagnostics produced by
the AHR-ES optimiser.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the fitted
parameter values stored in an `ahres.results` object.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) is currently a
placeholder method for `ahres.results` objects. No graphical summary is
produced at present.

## Usage

``` r
# S3 method for class 'ahres.results'
print(x, short = FALSE, ...)

# S3 method for class 'ahres.results'
coef(object, ...)

# S3 method for class 'ahres.results'
plot(x, ...)

# S3 method for class 'ahres.results'
summary(object, ...)

# S3 method for class 'summary.ahres.results'
print(x, ...)
```

## Format

A list, usually with the following components:

- par:

  Best parameter values found by the AHR-ES algorithm.

- value:

  Aggregated objective-function value at `par`.

- counts:

  Named integer vector giving the number of function and gradient
  evaluations recorded by the optimiser.

- convergence:

  Integer convergence code. `0` indicates that the stopping criterion
  was reached before the generation limit; `1` indicates that the
  maximum number of generations or function evaluations was reached.

- message:

  Character string describing the stopping condition.

- generations:

  Number of generations completed.

- partial:

  Optional vector of partial objective-function values, present when the
  objective function returns multiple components.

- trace:

  Optional tracing information collected during the run.

- active:

  Information on which parameters were active in the optimisation, as
  stored by
  [`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md).

## Arguments

- x:

  An `ahres.results` object, or a `summary.ahres.results` object for
  `print.summary.ahres.results()`.

- short:

  Logical. If `TRUE`, print a compact summary without listing the full
  parameter vector.

- ...:

  Additional arguments passed to the underlying method.

- object:

  An `ahres.results` object.

## Value

`coef.ahres.results()` returns the fitted parameter vector.

`plot.ahres.results()` returns `NULL`, invisibly.

`summary.ahres.results()` returns an object of class
`summary.ahres.results`.

## Details

These objects are specific to the AHR-ES optimisation workflow. They
contain the final parameter vector, the aggregated objective value,
convergence and iteration information, and, when relevant, partial
objective values and tracing information.

The main accessors are [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`coef()`](https://rdrr.io/r/stats/coef.html). A
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method is
currently exported as a placeholder but does not yet produce a graphical
summary.

## See also

[`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html)
