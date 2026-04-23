# Create small demonstration calibration projects

Create self-contained example projects that can be used to learn the
`calibrar` workflow, test installation, or prototype calibration code.
Each demo writes a small set of files to disk and returns the
information needed to read the setup table, load the observed data,
define the model interface, and run a calibration.

## Usage

``` r
calibrar_demo(path = NULL, model = NULL, ...)

# S3 method for class 'calibrar.demo'
print(x, ...)
```

## Arguments

- path:

  Directory where the demo project should be created. If `NULL`, the
  current working directory is used.

- model:

  Name of the demo model to generate. Supported values currently include
  `"PoissonMixedModel"`, `"PredatorPrey"`, `"SIR"`, and
  `"IBMLotkaVolterra"`.

- ...:

  Further arguments passed to
  [`print()`](https://rdrr.io/r/base/print.html).

- x:

  A `calibrar.demo` object returned by `calibrar_demo()`.

## Value

An object of class `calibrar.demo`, which also inherits from
[`calibrar.results()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md).
The returned list contains the demo project path, the true parameter
values used to generate the synthetic data, the path to the calibration
setup file, suggested starting values and bounds, phase information,
model-specific constants, and a precomputed objective value at the true
parameter set.

## Details

The generated project is meant to be a compact, runnable example of the
`calibrar` workflow rather than a realistic scientific application. The
returned object typically contains enough information to proceed through
the standard helper functions:

- [`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md)
  to read the setup table,

- [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md)
  to load the observed data,

- [`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
  to build the objective function, and

- [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
  or
  [`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
  to run the optimisation.

The exact files created depend on `model`, but they usually include a
calibration settings table, one or more observation files, and any
required forcing or master data files.

The currently available demos are:

- `"PoissonMixedModel"`:

  A Poisson autoregressive mixed model for the dynamics of a population
  observed at several sites over time.

- `"PredatorPrey"`:

  A deterministic Lotka–Volterra predator–prey system solved
  numerically.

- `"SIR"`:

  A simple susceptible–infected–recovered epidemic model solved
  numerically.

- `"IBMLotkaVolterra"`:

  A stochastic individual-based predator–prey example.

## See also

[`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md),
[`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md),
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md),
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md),
[`summary.calibrar.results()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md),
[`print.calibrar.results()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)

## Examples

``` r
tmp = tempdir()
demo = calibrar_demo(path = tmp, model = "PoissonMixedModel", T = 8, L = 3)
#> Creating observed data list for calibration...
#> Loaded observed data for variables: 'site_1', 'site_2', 'site_3'.
print(demo)
#> Calibration demo using model 'PoissonMixedModel'.
#> Target value: -6988.968 
#> Parameters:
#>        alpha         beta       gamma1       gamma2       gamma3       gamma4 
#>  0.400000000 -0.400000000 -0.011421355  0.100721594  0.217153872 -0.138190768 
#>       gamma5       gamma6       gamma7           sd      mu_ini1      mu_ini2 
#> -0.256919871  0.009345234 -0.047141311  0.200000000  4.330733340  4.624972813 
#>      mu_ini3 
#>  4.204692619 
#> Lower bounds:
#>   alpha    beta  gamma1  gamma2  gamma3  gamma4  gamma5  gamma6  gamma7      sd 
#>     0.0    -0.5    -1.0    -1.0    -1.0    -1.0    -1.0    -1.0    -1.0     0.0 
#> mu_ini1 mu_ini2 mu_ini3 
#>     3.7     3.7     3.7 
#> Upper bounds:
#>   alpha    beta  gamma1  gamma2  gamma3  gamma4  gamma5  gamma6  gamma7      sd 
#>     1.0     0.5     1.0     1.0     1.0     1.0     1.0     1.0     1.0     1.0 
#> mu_ini1 mu_ini2 mu_ini3 
#>     5.0     5.0     5.0 
#> Phases:
#>   alpha    beta  gamma1  gamma2  gamma3  gamma4  gamma5  gamma6  gamma7      sd 
#>       1       1       1       1       1       1       1       1       1      NA 
#> mu_ini1 mu_ini2 mu_ini3 
#>      NA      NA      NA 
#> Guess (starting point):
#>   alpha    beta  gamma1  gamma2  gamma3  gamma4  gamma5  gamma6  gamma7      sd 
#>   0.200   0.100   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.200 
#> mu_ini1 mu_ini2 mu_ini3 
#>   4.431   4.554   4.357 
#> 
#> Files location: C:\Users\ricar\AppData\Local\Temp\RtmpMT8afP/PoissonDemo

setup = calibration_setup(demo$setup)
observed = calibration_data(setup, path = demo$path)
#> Creating observed data list for calibration...
#> Loaded observed data for variables: 'site_1', 'site_2', 'site_3'.
```
