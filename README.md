# calibrar

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/calibrar)](http://cran.r-project.org/package=calibrar)
[![Github Issues](http://githubbadges.herokuapp.com/roliveros-ramos/calibrar/issues.svg?style=flat-square)](https://github.com/roliveros-ramos/calibrar/issues)

**Automated calibration for complex (ecological) models in R.**  
  The calibrar package allows the calibration of complex models, 
  including stochastic ones. It is a generic tool that can be used for 
  any type of model, especially those with non-differentiable objective functions. 
  It supports multiple phase calibrations and constrained optimization. 
  It implements maximum likelihood estimation methods and automated construction 
  of the objective function from simulated model outputs. 
  See <http://roliveros-ramos.github.io/calibrar> for more details.

## Installation

Get the released version from CRAN:

```R
install.packages("calibrar")
```

Or the development version from github:

```R
# install.packages("devtools")
devtools::install_github("roliveros-ramos/calibrar")
```
