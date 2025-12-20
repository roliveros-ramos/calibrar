# Calculate a discretization of the 2D Gaussian Kernel

Calculate a discretization of the 2D Gaussian Kernel

## Usage

``` r
gaussian_kernel(par, lower, upper, n = 10, checkSymmetry = TRUE, ...)
```

## Arguments

- par:

  A list, including the mean and covariance matrix.

- lower:

  A vector, indicating the lower bound for the calculation.

- upper:

  A vector, indicating the upper bound for the calculation.

- n:

  The number of cells for each dimension, can be one or two numbers.

- checkSymmetry:

  TRUE by default, checks if the covariance matrix is symmetric.

- ...:

  Additional arguments, currently not used.

## Value

A list, with 'x', 'y' and 'z' components.
