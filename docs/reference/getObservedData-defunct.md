# Get observed data for the calibration of a model

Create a list with the observed data with the information provided by
its main argument.

## Arguments

- info:

  A data.frame with the information about the calibration, normally
  created with the
  [`getCalibrationInfo`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  function. See details.

- path:

  Path to the directory to look up for the data.

- data.folder:

  folder in the path containing the data.

- ...:

  Additional arguments to `read.csv` function to read the data files.

## Value

A list with the observed data needed for a calibration, to be used in
combination with the
[`createObjectiveFunction`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md).

## See also

[`calibrar-defunct`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)

## Author

Ricardo Oliveros-Ramos
