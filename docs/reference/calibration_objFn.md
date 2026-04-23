# Create a calibration objective function from a model, setup table, and data

Build an objective function that runs a model, compares its outputs with
the observed data prepared by
[`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md),
and returns either a vector of partial objective values or a single
aggregated objective value.

## Usage

``` r
calibration_objFn(model, setup, observed, aggFn = NULL, aggregate = FALSE, ...)
```

## Arguments

- model:

  Function that runs the model for a candidate parameter set. Its first
  argument must be `par`. The function must return a named list of
  simulated outputs whose names match `setup$variable` for all active
  calibration rows.

- setup:

  A calibration setup table, usually created with
  [`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md).

- observed:

  Observed data as returned by
  [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md).
  This is typically a named list whose elements correspond to
  `setup$variable`.

- aggFn:

  Function used to aggregate partial objective values when
  `aggregate = TRUE`. The default is an internal weighted sum. The
  function must accept arguments `x` (partial objective values) and `w`
  (weights).

- aggregate:

  Logical. If `FALSE` (the default), the returned function produces one
  objective value per active calibration row. If `TRUE`, those partial
  values are combined into a single scalar using `aggFn` and the weights
  stored in `setup`.

- ...:

  Additional arguments passed to `model` each time the returned
  objective function is evaluated.

## Value

A function of `par` suitable for optimisation.

When evaluated, the returned function:

1.  runs `model(par, ...)`;

2.  extracts the simulated outputs;

3.  compares them with `observed` according to `setup$type`;

4.  returns either a named numeric vector of partial objective values or
    a single aggregated numeric value.

The returned function has class `"objFn"` and carries metadata
attributes used internally by
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
including the number of objective components, their weights, and the
names of the active variables.

## Details

This is the main bridge between a model-specific simulation function and
the optimisation tools in `calibrar`. The resulting function can be
passed directly to
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
and related optimisers.

For each row where `setup$calibrate` is `TRUE`, the function looks for a
variable named `setup$variable` in both `observed` and the list returned
by `model`. The comparison itself is delegated to the objective function
named in `setup$type`, using the same convention as
[`objFn()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md).

This means that the structure of the objects returned by `model` must
already be compatible with the corresponding objective function. For
example:

- pointwise criteria such as `"norm2"`, `"lnorm2"`, or `"pois"` usually
  expect vectors, matrices, or arrays with shapes matching the observed
  data;

- `"multinom"` expects matrix inputs;

- penalty terms such as `"normp"` may ignore `observed`, but still
  require the model output to contain the corresponding named element.

If `model` returns `NULL`, the generated objective function returns
`NULL` after emitting a diagnostic message.

## Aggregated versus vector objectives

With `aggregate = FALSE`, the returned objective is multi-component: it
produces one value per active calibration row. This can be useful for
methods that exploit vector-valued objectives.

With `aggregate = TRUE`, the returned objective is scalar and can be
used with optimisers that require a single value. Aggregation is
performed with
`aggFn(x = partial_values, w = setup$weight[setup$calibrate])`.

## Expected model interface

In typical use, the model function has the form:

`function(par, ...) -> named list`

where `par` may be a numeric vector or a list, depending on how the
optimisation is configured. The returned list should expose one element
per calibration target, with names that match `setup$variable` exactly.

## Workflow

The usual sequence is:

1.  define a model-running function;

2.  read the setup table with
    [`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md);

3.  read observed data with
    [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md);

4.  create the objective function with `calibration_objFn()`;

5.  call
    [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
    with the resulting function.

## See also

[`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md),
[`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md),
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`objFn()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md)

## Examples

``` r
if (FALSE) { # \dontrun{
demo = calibrar_demo(model = "PredatorPrey")
setup = calibration_setup(file = demo$setup)
observed = calibration_data(setup = setup, path = demo$path)

run_model = calibrar:::.PredatorPreyModel

obj_vec = calibration_objFn(
  model = run_model,
  setup = setup,
  observed = observed,
  T = demo$T,
  aggregate = FALSE
)

obj_scalar = calibration_objFn(
  model = run_model,
  setup = setup,
  observed = observed,
  T = demo$T,
  aggregate = TRUE
)

obj_vec(demo$guess)
obj_scalar(demo$guess)
} # }
```
