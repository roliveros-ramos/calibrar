# Create an objective function to be used with optimization routines

Create a new function, to be used as the objective function in the
calibration, given a function to run the model within R, observed data
and information about the comparison with data.

## Usage

``` r
calibration_objFn(model, setup, observed, aggFn = NULL, aggregate = FALSE, ...)
```

## Arguments

- model:

  Function to run the model and produce a list of outputs.

- setup:

  A data.frame with the information about the calibration, normally
  created with the
  [`calibration_setup`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md)
  function. See details.

- observed:

  A list of the observed variables created with the function
  [`calibration_data`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md)

- aggFn:

  A function to aggregate `fn` to a scalar value if the returned value
  is a vector. Some optimization algorithm can explote the additional
  information provided by a vectorial output from `fn`

- aggregate:

  boolean, if TRUE, a scalar value is returned using the `aggFn`.

- ...:

  More arguments passed to the `model` function.

## Value

A function, integrating the simulation of the model and the comparison
with observed data.

## See also

[`calibration_data`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md),
[`calibration_setup`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md).

## Author

Ricardo Oliveros-Ramos
