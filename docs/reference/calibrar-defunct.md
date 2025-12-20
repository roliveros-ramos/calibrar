# Defunct functions in package calibrar.

The functions listed below are defunct. When possible, alternative
functions with similar functionality are also mentioned. Help pages for
defunct functions are available at `help("<function>-defunct")`.

## Usage

``` r
getObservedData(info, path, data.folder = "data", ...)

getCalibrationInfo(
  path,
  file = "calibrationInfo.csv",
  stringsAsFactors = FALSE,
  ...
)

createObjectiveFunction(
  runModel,
  info,
  observed,
  aggFn = .weighted.sum,
  aggregate = FALSE,
  ...
)
```

## `getObservedData`

Deprecated in v0.3, defunct in v0.9. For `getObservedData`, use
[`calibration_data`](https://roliveros-ramos.github.io/calibrar/reference/calibration_data.md).

## `getCalibrationInfo`

Deprecated in v0.3, defunct in v0.9. For `getCalibrationInfo`, use
[`calibration_setup`](https://roliveros-ramos.github.io/calibrar/reference/calibration_setup.md).

## `createObjectiveFunction`

Deprecated in v0.3, defunct in v0.9. For `createObjectiveFunction`, use
[`calibration_objFn`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md).
