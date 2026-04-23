# Package index

## Main functions

Functions for general purpose optimization, with support for parallel
computation of derivatives.

- [`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
  : Unified optimisation interface with structured parameters and
  parallel numerical gradients
- [`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)
  : General-purpose optimisation with heuristic and derivative-free
  methods
- [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
  : Sequential parameter estimation for the calibration of complex
  models
- [`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md)
  : Adaptive Hierarchical Recombination Evolutionary Strategy
- [`gradient()`](https://roliveros-ramos.github.io/calibrar/reference/gradient.md)
  : Numerical gradient approximation
- [`calibrar-package`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-package.md)
  [`calibrar`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-package.md)
  : Automated Calibration for Complex Models

## Setting up a calibration

Functions to easily setup a new calibration for a complex model, from an
R function running the model.

- [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md)
  : Read observed data defined in a calibration setup table
- [`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
  : Create a calibration objective function from a model, setup table,
  and data
- [`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md)
  : Read and validate a calibration setup table
- [`objFn()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md)
  [`fitness()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md)
  : Objective function between observed and simulated data
- [`calibrar_demo()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar_demo.md)
  [`print(`*`<calibrar.demo>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar_demo.md)
  : Create small demonstration calibration projects

## Parameter modelling

Functions to easily parametrize a model.

- [`gaussian_kernel()`](https://roliveros-ramos.github.io/calibrar/reference/gaussian_kernel.md)
  : Discretise a bivariate Gaussian kernel on a grid
- [`spline_par()`](https://roliveros-ramos.github.io/calibrar/reference/spline_par.md)
  [`plot(`*`<spline.par>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/spline_par.md)
  : Interpolate time-varying parameters with a spline

## Test optimization functions

- [`sphereN()`](https://roliveros-ramos.github.io/calibrar/reference/sphereN.md)
  : Noisy sphere test function

## Auxiliar functions

- [`.get_command_argument()`](https://roliveros-ramos.github.io/calibrar/reference/dot-get_command_argument.md)
  : Extract a named command-line argument
- [`.read_configuration()`](https://roliveros-ramos.github.io/calibrar/reference/dot-read_configuration.md)
  : Read a simple key-value configuration file

## Defunct functions

- [`getObservedData()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  [`getCalibrationInfo()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  [`createObjectiveFunction()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  :

  Defunct functions in package calibrar.

## Result objects

- [`print(`*`<calibrar.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)
  [`coef(`*`<calibrar.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)
  [`predict(`*`<calibrar.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)
  [`plot(`*`<calibrar.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)
  [`summary(`*`<calibrar.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)
  [`print(`*`<summary.calibrar.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)
  [`dim(`*`<calibrar.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/calibrar.results.md)
  : Calibration result objects
- [`print(`*`<ahres.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/ahres.results.md)
  [`coef(`*`<ahres.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/ahres.results.md)
  [`plot(`*`<ahres.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/ahres.results.md)
  [`summary(`*`<ahres.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/ahres.results.md)
  [`print(`*`<summary.ahres.results>`*`)`](https://roliveros-ramos.github.io/calibrar/reference/ahres.results.md)
  : Evolutionary-strategy result objects
