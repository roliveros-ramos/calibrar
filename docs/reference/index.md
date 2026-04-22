# Package index

## Main functions

Functions for general purpose optimization, with support for parallel
computation of derivatives.

- [`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
  : Unified optimisation interface with structured parameters and
  parallel numerical gradients
- [`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)
  : General-purpose optimization using heuristic algorithms
- [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
  : Sequential parameter estimation for the calibration of complex
  models
- [`ahres()`](https://roliveros-ramos.github.io/calibrar/reference/ahres.md)
  : Adaptative Hierarchical Recombination Evolutionary Strategy (AHR-ES)
  for derivative-free and black-box optimization
- [`gradient()`](https://roliveros-ramos.github.io/calibrar/reference/gradient.md)
  : Numerical computation of the gradient, with parallel capabilities
- [`calibrar-package`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-package.md)
  [`calibrar`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-package.md)
  : Automated Calibration for Complex Models

## Setting up a calibration

Functions to easily setup a new calibration for a complex model, from an
R function running the model.

- [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md)
  : Get observed data for the calibration of a model

- [`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
  : Create an objective function to be used with optimization routines

- [`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md)
  :

  Get information to run a calibration using the `calibrar` package.

- [`objFn()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md)
  [`fitness()`](https://roliveros-ramos.github.io/calibrar/reference/objFn.md)
  : Objective function between observed and simulated data

- [`calibrar_demo()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar_demo.md)
  : Demos for the calibrar package

## Parameter modelling

Functions to easily parametrize a model.

- [`gaussian_kernel()`](https://roliveros-ramos.github.io/calibrar/reference/gaussian_kernel.md)
  : Calculate a discretization of the 2D Gaussian Kernel
- [`spline_par()`](https://roliveros-ramos.github.io/calibrar/reference/spline_par.md)
  : Predict time-varying parameters using splines.

## Test optimization functions

- [`sphereN()`](https://roliveros-ramos.github.io/calibrar/reference/sphereN.md)
  : Sphere function with random noise

## Auxiliar functions

- [`.get_command_argument()`](https://roliveros-ramos.github.io/calibrar/reference/dot-get_command_argument.md)
  : Get an specific argument from the command line
- [`.read_configuration()`](https://roliveros-ramos.github.io/calibrar/reference/dot-read_configuration.md)
  : Read a configuration file.

## Defunct functions

- [`getObservedData()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  [`getCalibrationInfo()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  [`createObjectiveFunction()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  :

  Defunct functions in package calibrar.
