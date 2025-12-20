# Create an objective function to be used with optimization routines

Create a new function, to be used as the objective function in the
calibration, given a function to run the model within R, observed data
and information about the comparison with data.

## Arguments

- runModel:

  Function to run the model and produce a list of outputs.

- info:

  A data.frame with the information about the calibration, normally
  created with the
  [`getCalibrationInfo`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  function. See details.

- observed:

  A list of the observed variables created with the function
  [`getObservedData`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)

- aggFn:

  A function to aggregate `fn` to a scalar value if the returned value
  is a vector. Some optimization algorithm can explote the additional
  information provided by a vectorial output from `fn`

- aggregate:

  boolean, if TRUE, a scalar value is returned using the `aggFn`.

- ...:

  More arguments passed to the `runModel` function.

## Value

A function, integrating the simulation of the model and the comparison
with observed data.

## See also

[`calibrar-defunct`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)

## Author

Ricardo Oliveros-Ramos
