# Read a simple key-value configuration file

Read a plain-text configuration file in which each non-empty,
non-comment line defines a key and one or more values separated by a
recognised delimiter. Parsed values are returned as a named list.

## Usage

``` r
.read_configuration(
  file,
  recursive = TRUE,
  keep.names = TRUE,
  conf.key = NULL,
  ...
)
```

## Arguments

- file:

  Path to the configuration file.

- recursive:

  Logical. If `TRUE`, keys matching `conf.key` are treated as references
  to additional configuration files and are read recursively.

- keep.names:

  Logical. If `TRUE`, keys are kept as written. If `FALSE`, they are
  converted to lower case.

- conf.key:

  Optional pattern identifying keys that should be interpreted as
  additional configuration files.

- ...:

  Reserved for future extensions.

## Value

A named list containing the parsed configuration values. Scalars are
returned as length-one vectors where appropriate, and simple values are
converted with
[`type.convert()`](https://rdrr.io/r/utils/type.convert.html). If `file`
does not exist, `NULL` is returned with a warning.

## Details

The reader accepts simple lines of the form `key separator value`, where
the separator may vary across lines. The current parser guesses the
separator among common choices such as `=`, `;`, `,`, `:`, and tab.

Lines beginning with punctuation after trimming are ignored, and
trailing comments are removed. When `conf.key` is supplied, matching
entries are treated as file references. Their paths are resolved
relative to the parent configuration file, and the referenced files can
be merged recursively.

This helper is intentionally lightweight. It is useful for small
calibration or script configuration files, but it is not intended to
replace richer configuration formats such as YAML.

## See also

[`.get_command_argument()`](https://roliveros-ramos.github.io/calibrar/reference/dot-get_command_argument.md)

## Examples

``` r
tf = tempfile(fileext = ".txt")
writeLines(c(
  "ncores = 4",
  "verbose = TRUE",
  "method = AHR-ES"
), tf)

cfg = .read_configuration(tf)
str(cfg)
#> List of 3
#>  $ ncores : int 4
#>  $ verbose: logi TRUE
#>  $ method : chr "AHR-ES"
```
