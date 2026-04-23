# Extract a named command-line argument

Parse command-line arguments and return the value associated with a
named option. This helper is primarily intended for scripting workflows
that call R from the shell and pass arguments such as
`--config=file.txt`.

## Usage

``` r
.get_command_argument(
  x,
  argument,
  prefix = "--",
  default = FALSE,
  verbose = FALSE
)
```

## Arguments

- x:

  Character vector of command-line arguments, typically obtained from
  [`commandArgs()`](https://rdrr.io/r/base/commandArgs.html).

- argument:

  Name of the argument to retrieve, without the prefix.

- prefix:

  Prefix used for arguments of interest. The default is `"--"`.

- default:

  Value returned when the argument is not found. The default is `FALSE`.

- verbose:

  Logical. If `TRUE`, a warning is issued when the requested argument is
  not found.

## Value

The parsed argument value. If the argument is present without an
explicit `=value` part, `TRUE` is returned. If the argument is absent,
`default` is returned. Parsed values are passed through an internal type
guesser, so logical and numeric values may be returned with their
corresponding R types.

## Details

The function looks for entries matching either `prefix + argument` or
`prefix + argument = value`. For example, with `prefix = "--"` and
`argument = "config"`, both `--config` and `--config=file.txt` are
recognised.

This function is exported because it can be useful in external packages
or scripts that build command-line interfaces around `calibrar`
workflows.

## See also

[`.read_configuration()`](https://roliveros-ramos.github.io/calibrar/reference/dot-read_configuration.md)

## Examples

``` r
.get_command_argument(commandArgs(), "interactive")
#> [1] FALSE
.get_command_argument(commandArgs(), "RStudio")
#> [1] FALSE
.get_command_argument(commandArgs(), "RStudio", prefix = "")
#> [1] FALSE
.get_command_argument(commandArgs(), "vanilla")
#> [1] FALSE
.get_command_argument("--control.file=baz.txt", "control.file")
#> [1] "baz.txt"
```
