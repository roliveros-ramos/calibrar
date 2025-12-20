# Changelog

## calibrar 0.9

- new
  [`optim2()`](https://roliveros-ramos.github.io/calibrar/reference/optim2.md)
  equivalent to [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  but with parallel computation of numerical gradients.

- new
  [`optimh()`](https://roliveros-ramos.github.io/calibrar/reference/optimh.md)
  wrapping heuristic methods with the same syntax of
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- the
  [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
  function implements the restart functionality for the `Rvmmin` method
  too, useful for the optimization of deterministic functions with long
  runtime.

- Improved methods for visualization of results.

- All optimization methods available in
  [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
  can use functions reading and writing from the disk.

- Function[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
  can use a different method for each estimation phase.

- [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
  is a generic now.

- Automatic stopping criteria for the AHR-ES method:

  ``` R
      - 0: maxit/maxgen only
      - 1: 1 OR max step reduction
      - 2: relative tolerance on value (smoothing for AHR-ES)
      - 3: maximum number of generations without improvement of `reltol`.
  ```

- Automatic testing using `testthat` package.

- Automatic support to optimize functions produced with the `TMB`
  package, via a method for
  [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md).

- [`getCalibrationInfo()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md),
  `createObjectiveFuction()` and
  [`getObservedData()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  are defunct now.

## calibrar 0.3

- new optimization methods available in
  [`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md):
  ‘LBFGSB3’, ‘hjn’, ‘CMA-ES’, ‘genSA’, ‘DE’, ‘soma’, ‘genoud’, ‘PSO’,
  ‘hybridPSO’, ‘mads’.
- fine control of numerical gradient computations, including
  parallelization.
- replicates argument for stochastic functions
- several minor bugs fixed
- [`getCalibrationInfo()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md),
  `createObjectiveFuction()` and
  [`getObservedData()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar-defunct.md)
  are deprecated and replaced by
  [`calibration_setup()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md),
  [`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
  and
  [`calibration_data()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md).
- [`spline_par()`](https://roliveros-ramos.github.io/calibrar/reference/spline_par.md)
  function to simplify the estimation of smooth time-varying parameters.

## calibrar 0.2

- par argument for the calibrate function can be a list
- optimization methods from `optimx`,
  [`stats::optim`](https://rdrr.io/r/stats/optim.html) and `cmaes` can
  be used
- several minor bugs fixed

## calibrar 0.1

- First release
