# Get observed data for the calibration of a model

Create a list with the observed data with the information provided by
its main argument.

## Usage

``` r
calibration_data(setup, path = ".", file = NULL, verbose = TRUE, ...)
```

## Arguments

- setup:

  A data.frame with the information about the calibration, normally
  created with the
  [`calibration_setup`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md)
  function. See details.

- path:

  Path to the directory to look up for the data. Paths in setup are
  considered relatives to this path.

- file:

  Optional file to save the created object (as an 'rds' file.)

- verbose:

  If TRUE, detailed messages of the process are printed.

- ...:

  Additional arguments to `read.csv` function to read the data files.

## Value

A list with the observed data needed for a calibration, to be used in
combination with the
[`calibration_objFn`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md).

## See also

[`calibration_objFn`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
[`calibration_setup`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md).

## Author

Ricardo Oliveros-Ramos
