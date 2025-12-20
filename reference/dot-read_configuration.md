# Read a configuration file.

File is expected to have lines of the form 'key SEP value' where key is
the name of the parameter, SEP a separator (can be '=' ',', ';') and
value the value of the parameter itself. The SEP for each line is
determined and parameters values are returned as a list.

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

  File to read the configuration

- recursive:

  Should 'conf.key' keys be read as additional configuration files?
  Default is TRUE.

- keep.names:

  Should names be kept as they are? By default, are converted to lower
  case.

- conf.key:

  String indicating the leading key to find an additional configuration
  file.

- ...:

  Additional arguments, not currently in use.
