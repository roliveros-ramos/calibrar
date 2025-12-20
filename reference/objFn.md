# Calcuted error measure between observed and simulated data

Calcuted error measure between observed and simulated data

## Usage

``` r
objFn(obs, sim, FUN, ...)

fitness(obs, sim, FUN, ...)
```

## Arguments

- obs:

  observed data as expected by FUN.

- sim:

  simulated data matching 'obs'

- FUN:

  the error function. Current accepted values area: 'norm2', 'lnorm2',
  'lnorm3', 'multinomial', 'pois', 'penalty0', 'penalty1', 'penalty2'
  and 'normp'.

- ...:

  Additional arguments to FUN

## Value

the value of FUN(obs, sim, ...)
