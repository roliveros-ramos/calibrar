# Get information to run a calibration using the `calibrar` package.

A wrapper for `read.csv` checking column names and data types for the
table with the calibration information.

## Arguments

- path:

  The path to look for the file.

- file:

  The file with the calibration information, see details.

- stringsAsFactors:

  To be passed to `read.csv`.

- ...:

  Additional arguments to `read.csv` function.

## Value

A data.frame with the information for the calibration of a model, to be
used with the
[`createObjectiveFunction`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
and
[`getObservedData`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md).

## See also

[`calibrar-defunct`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)

## Author

Ricardo Oliveros-Ramos
