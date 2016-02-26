# calibrar

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/calibrar)](http://cran.r-project.org/package=calibrar)
[![Github Issues](http://githubbadges.herokuapp.com/roliveros-ramos/calibrar/issues.svg?style=flat-square)](https://github.com/roliveros-ramos/calibrar/issues)
[![](http://cranlogs.r-pkg.org/badges/calibrar)](http://cran.rstudio.com/web/packages/calibrar/index.html)

**Automated parameter estimation for complex (ecological) models**  
  This package allows the parameter estimation or calibration of complex models, 
  including stochastic ones. It is a generic tool that can be used for fitting 
  any type of models, especially those with non-differentiable objective functions. 
  It supports multiple phases and constrained optimization. 
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
