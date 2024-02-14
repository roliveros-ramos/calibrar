# calibrar 0.9
* new `optim2()` equivalent to `stats::optim()` but with parallel computation of numerical gradients.
* new `optimh()` wrapping heuristic methods with the same syntax of `stats::optim()`.
* the `calibrate()` function implements the restart functionality for the `Rvmmin` method too, useful for the optimization of deterministic functions with long runtime.
* Improved methods for visualization of results.
* All optimization methods available in `calibrate()` can use functions reading and writing from the disk.
* Function`calibrate()` can use a different method for each estimation phase.
* `calibrate()` is a generic now.
* Automatic stopping criteria for the AHR-ES method:

          - 0: maxit/maxgen only
          - 1: 1 OR max step reduction
          - 2: relative tolerance on value (smoothing for AHR-ES)
          - 3: maximum number of generations without improvement of `reltol`.
* Automatic testing using `testthat` package.
* Automatic support to optimize functions produced with the `TMB` package, via a method for `calibrate()`.
* `getCalibrationInfo()`, `createObjectiveFuction()` and `getObservedData()` are defunct now.

# calibrar 0.3
* new optimization methods available in `calibrate()`: 'LBFGSB3', 'hjn', 'CMA-ES', 'genSA', 'DE', 'soma', 'genoud', 'PSO', 'hybridPSO', 'mads'.
* fine control of numerical gradient computations, including parallelization.
* replicates argument for stochastic functions 
* several minor bugs fixed
* `getCalibrationInfo()`, `createObjectiveFuction()` and `getObservedData()` are deprecated and replaced by `calibration_setup()`, `calibration_objFn()` and `calibration_data()`.
* `spline_par()` function to simplify the estimation of smooth time-varying parameters.


# calibrar 0.2
* par argument for the calibrate function can be a list
* optimization methods from `optimx`, `stats::optim` and `cmaes` can be used
* several minor bugs fixed

# calibrar 0.1
* First release

          
