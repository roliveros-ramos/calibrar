
<!-- README.md is generated from README.Rmd. Please edit that file -->

# calibrar <a href="https://roliveros-ramos.github.io/calibrar/"><img src="man/figures/logo_small.png" align="right" height="124" /></a>

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/calibrar)](https://CRAN.R-project.org/package=calibrar)
![GitHub R package
version](https://img.shields.io/github/r-package/v/roliveros-ramos/calibrar?label=GitHub)
[![R-CMD-check](https://github.com/roliveros-ramos/calibrar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/roliveros-ramos/calibrar/actions/workflows/R-CMD-check.yaml)
[![GitHub
issues](https://img.shields.io/github/issues/roliveros-ramos/calibrar)](https://github.com/roliveros-ramos/calibrar/issues)
[![OpenSSF Best
Practices](https://www.bestpractices.dev/projects/2132/badge)](https://www.bestpractices.dev/projects/2132)
[![](http://cranlogs.r-pkg.org/badges/calibrar)](https://CRAN.R-project.org/package=calibrar)
[![](http://cranlogs.r-pkg.org/badges/grand-total/calibrar)](https://CRAN.R-project.org/package=calibrar)
[![codecov](https://codecov.io/gh/roliveros-ramos/calibrar/graph/badge.svg?token=HELOL3WS4G)](https://app.codecov.io/gh/roliveros-ramos/calibrar)
<!-- badges: end -->

### Automated Parameter Estimation for Complex Models

This R package, calibrar, has been designed for the parameter estimation
(or calibration) of a wide range of ecological models, including complex
and stochastic models. The package combines various optimization
functionalities in a single interface, enabling the implementation of
the latest advancements in complex model calibration. The package
provides support for multiple phases and box constrained optimisation
with the possibility of using several algorithms available in R. In
particular, by using a “black-box” approach, the package allows the
calibration of models implemented in any programming language. It
provides a generic interface with models and allows the construction of
the objective function, within R, without requiring any changes in the
models’ code. Parallel support for computationally intensive models is
also provided, and can be used with high performance computing systems
in a simple manner, including the capability to restart an unfinished
optimisation for models with a long runtime.

It implements generic functions that can be used for fitting any type of
models, especially those with non-differentiable objective functions,
with the same syntax as base::optim. It supports multiple phases
estimation (sequential parameter masking), constrained optimization
(bounding box restrictions) and automatic parallel computation of
numerical gradients. Some common maximum likelihood estimation methods
and automated construction of the objective function from simulated
model outputs is provided.  
See <https://roliveros-ramos.github.io/calibrar/> for more details.

### Installation

``` r
# The easiest way to get calibrar is to install it from CRAN:
install.packages("calibrar")

# Alternatively, install the stable development version from OSMOSE drat repository:
install.packages("calibrar", repo="https://osmose-model.github.io/drat/")

# Or the development version from GitHub:
# install.packages("remotes")
remotes::install_github("roliveros-ramos/calibrar")
```

### Usage

For a quick introduction, check the worked the examples available from
the package:

``` r
library(calibrar)
vignette("calibrar")
```

For a more detailed explanation of the package philosophy, you can read
the article [calibrar: an R package for fitting complex ecological
models](https://doi.org/10.1111/2041-210X.14452).

### Contributions

If you find any bug, have questions about the documentation or requests
for enhancements, please [open an
issue](https://github.com/roliveros-ramos/calibrar/issues).

Contributions are accepted as pull requests. Please note that the
calibrar package is released with a [Contributor Code of
Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
By contributing to this project, you agree to abide by its terms.
