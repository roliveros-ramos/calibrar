# Get an specific argument from the command line

Get an specific argument from the command line

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

  The command line arguments, from `x = commandArgs()`

- argument:

  The name of the argument.

- prefix:

  The prefix to any argument of interest, the default is "–"

- default:

  Default value to return is argument is missing, default to FALSE.

- verbose:

  Boolean, if TRUE, shows a warning when the parameter is not found.

## Value

The value of the argument, assumed to be followed after '=' or, TRUE if
nothing but the argument was found. If the argument is not found, FALSE
is returned.

## Examples

``` r
.get_command_argument(commandArgs(), "interactive")
#> [1] FALSE
.get_command_argument(commandArgs(), "RStudio")
#> [1] FALSE
.get_command_argument(commandArgs(), "RStudio", prefix="")
#> [1] FALSE
.get_command_argument(commandArgs(), "vanilla")
#> [1] FALSE
.get_command_argument("--control.file=baz.txt", "control.file")
#> [1] "baz.txt"
```
