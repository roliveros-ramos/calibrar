# Get information to run a calibration using the `calibrar` package.

A wrapper for `read.csv` checking column names and data types for the
table with the calibration information.

## Usage

``` r
calibration_setup(file, control = list(), ...)
```

## Arguments

- file:

  The file with the calibration information, see details.

- control:

  Control arguments for generating the setup. See details.

- ...:

  Additional arguments to `read.csv` function.

## Value

A data.frame with the information for the calibration of a model, to be
used with the
[`calibration_objFn`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
and
[`calibration_data`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md).

## See also

[`calibration_objFn`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
[`calibration_data`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md).

## Author

Ricardo Oliveros-Ramos
