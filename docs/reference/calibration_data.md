# Read observed data defined in a calibration setup table

Build the `observed` object used by
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
from a setup table created with
[`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md).
Each active row in the setup table points to a file containing observed
data, and this function reads those files into a named list whose
elements correspond to the calibration variables.

## Usage

``` r
calibration_data(setup, path = ".", file = NULL, verbose = TRUE, ...)
```

## Arguments

- setup:

  A calibration setup table, usually created with
  [`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md).

- path:

  Base directory used to resolve relative paths from `setup$file`. Each
  file is looked up with
  [`base::file.path()`](https://rdrr.io/r/base/file.path.html) using
  `path` and the value stored in the corresponding `file` row.

- file:

  Optional path to save the resulting `observed` object with
  [`base::saveRDS()`](https://rdrr.io/r/base/readRDS.html). This is
  useful when observed data are expensive to read repeatedly.

- verbose:

  Logical. If `TRUE`, print progress messages describing which variables
  were loaded.

- ...:

  Additional arguments passed to
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html) when
  reading text files.

## Value

A named list with one element per row of `setup`. Element names are
taken from `setup$variable`.

Depending on the corresponding setup row, each element is typically one
of:

- a numeric vector, when `varid` selects a single column from a text
  file;

- a numeric matrix, when a full table is read from a text file;

- an arbitrary R object read from an `.rds` file;

- `NA`, when `use_data = FALSE` or `calibrate = FALSE`.

## Details

Rows with `use_data = FALSE` are kept in the output as `NA`. This makes
it possible to include penalty-only components in the same calibration
workflow: such rows still participate in
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
provided the model returns a simulated component with the same name.

Only rows satisfying both `setup$calibrate` and `setup$use_data` are
actually read from disk. All other rows are represented as `NA` in the
returned list.

File reading currently depends on file type:

- text files (`.csv`, `.txt`):

  Read with
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html). If
  `varid` is supplied, that column is extracted and returned as a
  numeric vector. Otherwise the full table is converted to a numeric
  matrix after dropping the first `col_skip` columns, when applicable.

- R objects (`.rds`):

  Read with [`base::readRDS()`](https://rdrr.io/r/base/readRDS.html) and
  returned as-is.

- netCDF (`.nc`):

  Not currently supported by the implementation.

The interpretation of `varid`, `col_skip`, and `nrows` is controlled by
the corresponding columns in `setup`, after defaults have been filled by
[`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md).

## Expected alignment with model outputs

The names of the returned list are used later by
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
to match observed and simulated quantities. For each active calibration
row, `setup$variable` must therefore match the name of the corresponding
element in the list returned by the model function.

## Saving and reusing observed data

When `file` is supplied, the created list is also saved as an `.rds`
file. This can be helpful in iterative calibration workflows where the
setup table is stable but reading the original data files repeatedly is
inconvenient.

## See also

[`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md),
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
[`objFn()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md)

## Examples

``` r
if (FALSE) { # \dontrun{
demo = calibrar_demo(model = "PredatorPrey")
setup = calibration_setup(file = demo$setup)
observed = calibration_data(setup = setup, path = demo$path)

names(observed)
str(observed)
} # }
```
