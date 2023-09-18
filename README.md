
<!-- README.md is generated from README.Rmd. Please edit that file -->

# calibrar <a href="https://roliveros-ramos.github.io/calibrar/"><img src="man/figures/logo_small.png" align="right" height="138" /></a>

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/calibrar)](http://cran.r-project.org/package=calibrar)
![GitHub R package
version](https://img.shields.io/github/r-package/v/roliveros-ramos/calibrar?label=GitHub)
[![GitHub
issues](https://img.shields.io/github/issues/roliveros-ramos/calibrar)](https://github.com/roliveros-ramos/calibrar/issues)
[![](http://cranlogs.r-pkg.org/badges/calibrar)](http://cran.rstudio.com/web/packages/calibrar/index.html)
[![OpenSSF Best
Practices](https://www.bestpractices.dev/projects/2132/badge)](https://www.bestpractices.dev/projects/2132)
<!-- badges: end -->

## Overview

This package allows the parameter estimation (i.e. calibration) of
complex models, including stochastic ones. It implements generic
functions that can be used for fitting any type of models, especially
those with non-differentiable objective functions, with the same syntax
as base::optim. It supports multiple phases estimation (sequential
parameter masking), constrained optimization (bounding box restrictions)
and automatic parallel computation of numerical gradients. Some common
maximum likelihood estimation methods and automated construction of the
objective function from simulated model outputs is provided.  
See <http://roliveros-ramos.github.io/calibrar> for more details.

## Installation

``` r
# The easiest way to get calibrar is to install it from CRAN:
install.packages("calibrar")

# Alternatively, install the stable development version from OSMOSE drat repository:
install.packages("calibrar", repo="https://osmose-model.github.io/drat/")

# Or the development version from GitHub:
# install.packages("pak")
remotes::install_github("roliveros-ramos/calibrar")
```

## Usage

For a quick introduction, check the worked the examples available from
the package:

``` r
library(calibrar)
vignette("calibrar")
```

For a more detailed explanation of the package philosophy, you can read
the pre-print [calibrar: an R package for fitting complex ecological
models](https://arxiv.org/abs/1603.03141).

## Contributions

If you find any bug, have questions about the documentation or requests
for enhancements, please [open an
issue](https://github.com/roliveros-ramos/calibrar/issues).

Contributions are accepted as pull requests. Please note that the
calibrar package is released with a [Contributor Code of
Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
By contributing to this project, you agree to abide by its terms.
