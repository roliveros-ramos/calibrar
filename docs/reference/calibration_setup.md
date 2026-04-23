# Read and validate a calibration setup table

Read a calibration settings table from a text file and coerce it to the
structure expected by
[`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md)
and
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md).
The setup table defines which model outputs are compared with
observations, which objective function is used for each output, whether
a row contributes to the calibration, and how each component is weighted
in the final objective.

## Usage

``` r
calibration_setup(file, control = list(), ...)
```

## Arguments

- file:

  Path to the calibration setup file, typically a `.csv` file.

- control:

  Named list of defaults used when optional columns are missing from the
  setup table. Currently recognised entries are:

  `col_skip`

  :   Default number of leading columns to discard when a text file is
      read as a full table by
      [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md).
      Defaults to `2`.

  `nrows`

  :   Default number of rows to read from each data file. Defaults to
      `NA`, which is interpreted downstream as “read all rows”.

- ...:

  Additional arguments passed to
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html) when
  reading `file`.

## Value

A `data.frame` describing the calibration targets. The returned table is
suitable for direct use in
[`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md)
and
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md).
Standard columns are coerced to the expected types and missing optional
columns are added where needed.

## Details

The function is a validated wrapper around
[`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html). It
checks required column names, coerces standard column types, fills a few
optional columns from `control`, and derives `weight` from `cv` when
needed.

The setup table may contain the following columns:

- `variable`:

  Name of the calibration component. This name must match the
  corresponding element returned by the model function supplied to
  [`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md).

- `type`:

  Name of the objective function used for that component. In current
  workflows this is usually the name of a function accepted by
  [`objFn()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md),
  such as `"norm2"`, `"lnorm2"`, `"pois"`, `"multinom"`, `"normp"`, or
  `"penalty"`.

- `calibrate`:

  Logical flag indicating whether the component is active in the
  calibration. Only rows with `calibrate = TRUE` are included in the
  objective function returned by
  [`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md).

- `weight`:

  Positive weight used when partial objective values are aggregated. If
  `weight` is missing, or contains `NA`, it is derived from `cv` as
  `1 / (2 * cv^2)` for the corresponding rows.

- `cv`:

  Coefficient of variation used to derive `weight` when an explicit
  weight is not supplied. At least one of `weight` or `cv` must be
  available in the table.

- `use_data`:

  Logical flag indicating whether an observed dataset must be read for
  that row. This is typically `FALSE` for pure penalties such as
  `"normp"`, where the objective depends only on simulated values.

- `file`:

  Path to the observed data file for that component. Paths are
  interpreted relative to the `path` argument of
  [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md).
  This column is still required when `use_data = FALSE`, although the
  file will not be read in that case.

- `varid`:

  Optional identifier used by
  [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md)
  when extracting a single column from a text file. If omitted, a full
  numeric table is read instead.

- `col_skip`:

  Optional number of leading columns to discard when a text file is read
  as a full table. Missing values are replaced by `control$col_skip`.

- `nrows`:

  Optional number of rows to read from each observed-data file. Missing
  values are replaced by `control$nrows`.

Required columns are `variable`, `type`, `calibrate`, `use_data`, and
`file`. In addition, at least one of `weight` or `cv` must be supplied.

## Weight handling

Aggregation weights are stored in the returned table regardless of
whether the objective is later used in vector or scalar form.

If the table does not contain a `weight` column, weights are computed
from `cv` as `1 / (2 * cv^2)`. If `weight` exists but contains missing
values, those entries are filled using the same formula. All final
weights must be strictly positive.

## Relationship with the rest of the workflow

A typical file-based calibration workflow is:

1.  read the setup table with `calibration_setup()`;

2.  load observed data with
    [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md);

3.  build the objective function with
    [`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md);

4.  pass the resulting function to
    [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
    or another optimiser.

## See also

[`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md),
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`objFn()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md)

## Examples

``` r
if (FALSE) { # \dontrun{
demo = calibrar_demo(model = "PredatorPrey")
setup = calibration_setup(file = demo$setup)
setup
} # }
```
