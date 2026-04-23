# Calibration result objects

Objects of class `calibrar.results` are returned by
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md).
They store the final parameter values together with the objective value,
convergence information, the set of calibrated parameters, and tracing
information from the calibration run.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the fitted
parameter values stored in a `calibrar.results` object.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) is currently a
placeholder method for `calibrar.results` objects. No graphical summary
is produced at present.

[`dim()`](https://rdrr.io/r/base/dim.html) returns the number of fitted
parameters after unlisting the parameter structure stored in the result
object.

## Usage

``` r
# S3 method for class 'calibrar.results'
print(x, ...)

# S3 method for class 'calibrar.results'
coef(object, ...)

# S3 method for class 'calibrar.results'
predict(object, replicates = 1, na.rm = FALSE, ...)

# S3 method for class 'calibrar.results'
plot(x, ...)

# S3 method for class 'calibrar.results'
summary(object, ..., show_par = NULL, par.only = FALSE)

# S3 method for class 'summary.calibrar.results'
print(x, digits = 3, ...)

# S3 method for class 'calibrar.results'
dim(x)
```

## Format

A list, usually with the following components:

- par:

  Best parameter values found, returned with the same structure as the
  input `par` supplied to
  [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md).

- value:

  Objective-function value at `par`.

- counts:

  Named integer vector giving the number of function and gradient
  evaluations recorded for the last calibration phase.

- convergence:

  Integer convergence code returned by the optimiser used in the final
  phase.

- message:

  Character string with additional information returned by the
  optimiser.

- method:

  Name of the optimisation method used in the calibration.

- fn:

  Objective function used for the calibration.

- active:

  Logical vector indicating which parameters were active (estimated) in
  the full calibration problem.

- elapsed:

  Elapsed time for the complete calibration run.

- trace:

  A list containing tracing information from the last phase and
  additional run-level metadata. Depending on the control settings and
  the optimiser used, this may include per-phase outputs, generation
  traces, partial objective values, and other diagnostic information.

## Arguments

- x:

  A `calibrar.results` object, or a `summary.calibrar.results` object
  for `print.summary.calibrar.results()`.

- ...:

  Additional arguments passed to the underlying method.

- object:

  A `calibrar.results` object.

- replicates:

  Number of replicate evaluations of the fitted parameter set to
  perform. If `replicates = 1`, the underlying model function is
  evaluated once. If `replicates > 1`, repeated evaluations are combined
  by taking the mean across replicates for each returned output.

- na.rm:

  Logical. Passed to the internal averaging step used when
  `replicates > 1`.

- show_par:

  Optional character vector of parameter names or integer vector of
  positions indicating which parameters should be shown in the summary.
  By default, all fitted parameters are shown.

- par.only:

  Logical. If `TRUE`, return only the method name and selected parameter
  values. This is mainly useful when summarising several
  `calibrar.results` objects side by side.

- digits:

  Number of significant digits to use when printing summary tables or
  parameter values.

## Value

`coef.calibrar.results()` returns the fitted parameter values, with the
same structure as the `par` component of the object.

`predict.calibrar.results()` returns the model outputs produced by the
fitted parameter set. When `replicates > 1`, it returns the mean output
across replicates and stores the full replicate-level outputs in the
`.replicates` component of the returned object.

`plot.calibrar.results()` returns `NULL`, invisibly.

`summary.calibrar.results()` returns either:

- an object of class `summary.calibrar.results` when called on a single
  object; or

- a `data.frame` with one row per result object when several
  `calibrar.results` objects are supplied through `...`.

The print method formats these summaries for console display.

`dim.calibrar.results()` returns a single integer.

## Details

A `calibrar.results` object represents the final output of a complete
calibration workflow, potentially including multiple calibration phases.

The most useful accessors are:

- [`print()`](https://rdrr.io/r/base/print.html) for a compact console
  summary;

- [`summary()`](https://rdrr.io/r/base/summary.html) for a more explicit
  summary, including selected parameters;

- [`coef()`](https://rdrr.io/r/stats/coef.html) to extract the fitted
  parameter values;

- [`predict()`](https://rdrr.io/r/stats/predict.html) to evaluate the
  fitted parameter set with the underlying model function when this is
  available; and

- [`dim()`](https://rdrr.io/r/base/dim.html) to report the number of
  parameters in the fitted solution.

The [`predict()`](https://rdrr.io/r/stats/predict.html) method is
currently only available when the objective function was built with
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
because it relies on metadata attached to that function.

## See also

[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`predict()`](https://rdrr.io/r/stats/predict.html)
