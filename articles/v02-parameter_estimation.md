# Using the \`calibrate()\` function for parameter estimation

## Introduction

This vignette focus on the use of the
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
for parameter estimation. We suggest to see the vignette ‘Getting
started with the `calibrar` package’ before reading this one, specially
if you do not have previous experience doing optimization in R.

## Estimating parameters for a linear model

As a first example, we will estimate the parameters for a linear model
by manually performing the optimization, in opposition of the standard
method using the [`stats::lm()`](https://rdrr.io/r/stats/lm.html). The
objetive of this is to introduce the features of the
[`calibrate()`](https://roliveros-ramos.github.io/calibrar/reference/calibrate.md)
function with a simple and fast model. Let’s start by creating some
parameters for the linear model.

``` r
library(calibrar)
N = 7 # number of variables in the linear model
T = 100 # number of observations
sd = 0.25 # standard deviation of the gaussian noise
# observed data
x = matrix(rnorm(N*T, sd=sd), nrow=T, ncol=N)
# slopes for the linear model (real parameters)
slope = seq_len(N) 
# intercept for the linear model (real parameters)
intercept = pi
# real parameters
real = list(intercept=intercept, slope=slope)
real
#> $intercept
#> [1] 3.141593
#> 
#> $slope
#> [1] 1 2 3 4 5 6 7
```

Now, let’s create a function so simulate the linear model.

``` r
# function to simulate the linear model
linear = function(x, par) {
  stopifnot(length(x)==length(par$slope))
  out = sum(x*par$slope) + par$intercept
  return(out)
}
```

And, finally, the simulated data for the exercise:

``` r
# simulated data 
y = apply(x, 1, linear, par=real)
```

Of course, the solution can be found using the
[`lm()`](https://rdrr.io/r/stats/lm.html) function:

``` r
mod = lm(y ~ x)
mod
#> 
#> Call:
#> lm(formula = y ~ x)
#> 
#> Coefficients:
#> (Intercept)           x1           x2           x3           x4           x5  
#>       3.142        1.000        2.000        3.000        4.000        5.000  
#>          x6           x7  
#>       6.000        7.000
```

Now, in order to proceed to find the solution by an explicit numerical
optimization, we need to define the objective function to be minimized:

``` r
# objective function (residual squares sum)
obj = function(par, x, y) {
  y_sim = apply(x, 1, linear, par=par)
  out = sum((y_sim - y)^2)
  return(out)
}
```

So now we can proceed with the optimization:

``` r
# initial guess for optimization
start = list(intercept=0, slope=rep(0, N))
bfgs = calibrate(par=start, fn=obj, x=x, y=y)
#> Using optimization method 'Rvmmin'.
#> Elapsed time: 0.15s
#> Function value: 5.28185e-14
#> Parameter values: 3.14 1 2 3 4 5 6 7
#> 
#> Status: Rvmminu appears to have converged
# using coef to extract optimal parameters
coef(bfgs)
#> $intercept
#> [1] 3.141593
#> 
#> $slope
#> [1] 1 2 3 4 5 6 7
```

As expected, we were able to recover the real parameters of the model.
Now, let’s specify `lower` and `upper` bounds for the algorithms that
require them:

``` r
lower = relist(rep(-10, N+1), skeleton=start)
upper = relist(rep(+10, N+1), skeleton=start)
```

And repeat the exercise with several optimization algorithms:

``` r
set.seed(880820) # for reproducibility
cg = calibrate(par=start, fn=obj, x=x, y=y, method='CG')
#> Using optimization method 'CG'.
#> Elapsed time: 0.56s
#> Function value: 3.37729e-13
#> Parameter values: 3.14 1 2 3 4 5 6 7
#> 
#> Status: -
nm = calibrate(par=start, fn=obj, x=x, y=y, method='nmkb', lower=lower, upper=upper)
#> Using optimization method 'nmkb'.
#> Elapsed time: 0.32s
#> Function value: 5.28315e-07
#> Parameter values: 3.14 1 2 3 4 5 6 7
#> 
#> Status: -
ahres = calibrate(par=start, fn=obj, x=x, y=y, method='AHR-ES')
#> Using optimization method 'AHR-ES'.
#> Elapsed time: 1.98s
#> Function value: 2.66292e-18
#> Parameter values: 3.14 1 2 3 4 5 6 7
#> 
#> Status: Stopping criteria reached in 234 generations.
hjn = calibrate(par=start, fn=obj, x=x, y=y, method='hjn', lower=lower, upper=upper)
#> Using optimization method 'hjn'.
#> Elapsed time: 0.52s
#> Function value: 2.15391e-13
#> Parameter values: 3.14 1 2 3 4 5 6 7
#> 
#> Status: hjn appears to have converged
```

And compare the results:

``` r
summary(ahres, hjn, nm, bfgs, cg, par.only=TRUE)
#>       method intercept slope1 slope2 slope3 slope4 slope5 slope6 slope7
#> ahres AHR-ES      3.14      1      2      3      4      5      6      7
#> hjn      hjn      3.14      1      2      3      4      5      6      7
#> nm      nmkb      3.14      1      2      3      4      5      6      7
#> bfgs  Rvmmin      3.14      1      2      3      4      5      6      7
#> cg        CG      3.14      1      2      3      4      5      6      7
```

As we can see, for this simple example, all the algorithms used were
able to find the solution within a reasonable time, but some were faster
than others. In the next examples we will see this is not always the
case, and some algorithms can perform very differently or even fail for
a particular optimization problem.

## Fitting a biomass production model with harvest

As a second example, we will estimate the parameters for a difference
equation system, simulating the dynamics of the biomass B of a harvested
population:

B\_{t+1} = B_t + rB_t\left(1-\frac{B_t}{K}\right) - C_t, where B_t is
the biomass in time t, C_t is the catch during the interval \[t, t+1\[,
r is the intrinsic population growth rate and K is the carrying capacity
of the system.

We will define some values for all the parameters so we can perform an
optimization an try to recover them from the simulated data:

``` r
set.seed(880820)
T = 50
real = list(r=0.5, K=1000, B0=600)
catch = 0.25*(real$r*real$K)*runif(T, min=0.2, max=1.8)
```

As before, a function to simulate data from the parameters (the model)
will be needed. The requirement for this function is to have as first
argument the parameter vector (or list) `par`:

``` r
run_model = function(par, T, catch) {
  B = numeric(T+1)
  times = seq(0, T)
  B0 = par$B0
  r = par$r
  K = par$K
  B[1] = B0
  for(t in seq_len(T)) {
    b = B[t] + r*B[t]*(1-B[t]/K) - catch[t] # could be negative
    B[t+1] = max(b, 0.01*B[t]) # smooth aproximation to zero
  }
  out = list(biomass=B)
  return(out)
}
```

And now we can use the `run_model()` function and the assumed parameters
to simulate the model:

``` r
observed = run_model(par=real, T=T, catch=catch)
```

``` r
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(1,1,1,1))
plot(observed$biomass, type="l", lwd=2, ylab="biomass", xlab="", las=1, ylim=c(0, 1.2*max(observed$biomass)))
mtext("BIOMASS", 3, adj=0.01, line = 0, font=2)
plot(catch, type="h", lwd=2, ylab="catch", xlab="", las=1, ylim=c(0, 1.2*max(catch)))
mtext("CATCH", 3, adj=0.01, line = 0, font=2)
```

![Simulation of the logistic models with the assumed
parameters.](v02-parameter_estimation_files/figure-html/logistic4-1.png)

In order to carry out the optimization, we need the objective function
to be defined, in this case, using a simple residual squares sum:

``` r
objfn = function(par, T, catch, observed) {
  simulated = run_model(par=par, T=T, catch=catch)
  value = sum((observed$biomass-simulated$biomass)^2, na.rm=TRUE)
  return(value)
}
```

Finally, we need to define a starting point for the search,

``` r
start = list(r=0.1, K=1.5*max(observed$biomass), B0=observed$biomass[1])
```

and we are ready to try to estimate the parameters using several
algorithms:

``` r
set.seed(880820) # for reproducibility
opt0 = calibrate(par=start, fn = objfn, method='LBFGSB3', T=T, catch=catch, observed=observed)
#> Using optimization method 'LBFGSB3'.
#> Elapsed time: 0.02s
#> Function value: 150998
#> Parameter values: 0.402 1.18e+03 600
#> 
#> Status: CONVERGENCE: Parameters differences below xtol
opt1 = calibrate(par=start, fn = objfn, method='Rvmmin', T=T, catch=catch, observed=observed)
#> Using optimization method 'Rvmmin'.
#> Elapsed time: 0.01s
#> Function value: 1.38129e+07
#> Parameter values: 2.9 1.14e+03 482
#> 
#> Status: Rvmminu appears to have converged
opt2 = calibrate(par=start, fn = objfn, method='CG', T=T, catch=catch, observed=observed)
#> Using optimization method 'CG'.
#> Elapsed time: 0.10s
#> Function value: 150931
#> Parameter values: 0.402 1.18e+03 600
#> 
#> Status: -
opt3 = calibrate(par=start, fn = objfn, method='AHR-ES', T=T, catch=catch, observed=observed)
#> Using optimization method 'AHR-ES'.
#> Elapsed time: 2.13s
#> Function value: 6.77545e-17
#> Parameter values: 0.5 1e+03 600
#> 
#> Status: Stopping criteria reached in 1590 generations.
opt4 = calibrate(par=start, fn = objfn, method='CMA-ES', T=T, catch=catch, observed=observed)
#> Using optimization method 'CMA-ES'.
#> Elapsed time: 0.18s
#> Function value: 0.914125
#> Parameter values: 0.5 1e+03 600
#> 
#> Status: Covariance matrix 'C' is numerically not positive definite.
opt5 = calibrate(par=start, fn = objfn, method='hjn', T=T, catch=catch, observed=observed)
#> Using optimization method 'hjn'.
#> Elapsed time: 0.54s
#> Function value: 990.8
#> Parameter values: 0.509 988 602
#> 
#> Status: Function count limit exceeded.
opt6 = calibrate(par=start, fn = objfn, method='Nelder-Mead', T=T, catch=catch, observed=observed)
#> Using optimization method 'Nelder-Mead'.
#> Elapsed time: 0.02s
#> Function value: 0.0775804
#> Parameter values: 0.5 1e+03 600
#> 
#> Status: -
```

The function [`summary()`](https://rdrr.io/r/base/summary.html) can be
used to compare all the optimization results:

``` r
summary(opt0, opt1, opt2, opt3, opt4, opt5, opt6)
#>           method elapsed    value    fn  gr     r    K  B0
#> opt0     LBFGSB3  0.0161 1.51e+05    18  18 0.402 1181 600
#> opt1      Rvmmin  0.0139 1.38e+07   108   6 2.900 1145 482
#> opt2          CG  0.1035 1.51e+05   691 101 0.402 1181 600
#> opt3      AHR-ES  2.1296 6.78e-17 11130   0 0.500 1000 600
#> opt4      CMA-ES  0.1838 9.14e-01  1540  NA 0.500 1000 600
#> opt5         hjn  0.5382 9.91e+02  6000  NA 0.509  988 602
#> opt6 Nelder-Mead  0.0208 7.76e-02   221  NA 0.500 1000 600
```

For a better comparison, we can simulate the results for all the
parameters found:

``` r
sim0 = run_model(par=coef(opt0), T=T, catch=catch)
sim1 = run_model(par=coef(opt1), T=T, catch=catch)
sim2 = run_model(par=coef(opt2), T=T, catch=catch)
sim3 = run_model(par=coef(opt3), T=T, catch=catch)
sim4 = run_model(par=coef(opt4), T=T, catch=catch)
sim5 = run_model(par=coef(opt5), T=T, catch=catch)
```

And plot some of the best results obtained (L-BFGS-B 3.0, CG and
AHR-ES):

``` r
par(mar=c(3,4,1,1))
plot(observed$biomass, type="n", ylab="BIOMASS", xlab="", las=1, ylim=c(0, 1.2*max(observed$biomass)))
lines(sim0$biomass, col=1, lwd=2)
lines(sim2$biomass, col=2, lwd=2)
lines(sim3$biomass, col=3, lwd=2)
points(observed$biomass)
mtext(c('LBFGSB3', 'CG', 'AHR-ES'), 1, adj=0.05, col=1:3, line=-(4:2), font=2)
```

![Plot of best results after parameter optimisation.
](v02-parameter_estimation_files/figure-html/unnamed-chunk-5-1.png)

So far, we have used all algorithms with their default arguments. Most
algorithms provide control arguments that allow to improve its
performance for a particular problem. For example, looking back to the
results from the L-BFGS-B 3.0 method, we can see as status ‘Maximum
number of iterations reached’, meaning the algorithm did not converge
but stopped after the maximum number of 100 iterations allowed by
default. This can be modified here (and for several other methods) by
changing the `maxit` control argument:

``` r
optx = calibrate(par=start, fn = objfn, method='LBFGSB3', T=T, catch=catch, observed=observed, control=list(maxit=20000))
#> Using optimization method 'LBFGSB3'.
#> Elapsed time: 0.01s
#> Function value: 150998
#> Parameter values: 0.402 1.18e+03 600
#> 
#> Status: CONVERGENCE: Parameters differences below xtol
```

And we can see that with this setup, the algorithm converges to the
right solution.

We also have the possibility to set up a calibration in multiple phases,
meaning we will solve several sequential optimizations with a
progressively higher number of parameters. The purpose of this is to
improve the initial search point for a final optimization with all the
parameters active. This heuristic may help to achieve find the solution
for some problems. For example, the Rvmmin algorithm has a status
‘Rvmminu appears to have converged’ but was not able to converge to the
original parameter values. So, we will try again by fixing some of the
parameters (the initial biomass) and trying a two phases parameter
estimation:

``` r
calibrate(par=start, fn = objfn, method='Rvmmin', T=T, catch=catch, observed=observed, phases = c(1,1,2))
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 3 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 1 finished (0.02s)
#>  Function value: 1.8303e-08
#>  Parameter values: 0.5 1e+03
#> 
#> - Phase 2: 3 out of 3 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 2 finished (0.00s)
#>  Function value: 1.8303e-08
#>  Parameter values: 0.5 1e+03 600
#> 
#> Status: Rvmminu appears to have converged
#> Optimization using 'Rvmmin' algorithm.
#> Function value: 1.830304e-08 
#> Status: Rvmminu appears to have converged 
#> Parameters:
#>     r     K    B0 
#> 5e-01 1e+03 6e+02 
#> Computation:
#> function gradient 
#>       17        1
```

And with two phases, now we also find the solution using the ‘Rvmmin’
method.

## Fitting an autoregressive Poisson model

As a third example, we will estimate the parameters for a Poisson
Autoregressive Mixed model for the dynamics of a population in different
sites:

log(\mu\_{i, t+1}) = log(\mu\_{i, t}) + \alpha + \beta X\_{i, t} +
\gamma_t

where \mu\_{i, t} is the size of the population in site i at year t,
X\_{i, t} is the value of an environmental variable in site i at year t.
The parameters to estimate were \alpha, \beta, and \gamma_t, the random
effects for each year, \gamma_t \sim N(0,\sigma^2), and the initial
population at each site \mu\_{i, 0}. We assumed that the observations
N\_{i,t} follow a Poisson distribution with mean \mu\_{i, t}. We could
also create the data for this model using the function
[`calibrar_demo()`](https://roliveros-ramos.github.io/calibrar/reference/calibrar_demo.md),
with the additional arguments `L=5` (five sites) and `T=100` (one
hundred years):

``` r
path = NULL # NULL to use the current directory
ARPM = calibrar_demo(path=path, model="PoissonMixedModel", L=5, T=100) 
#> Creating observed data list for calibration...
#> Loaded observed data for variables: 'site_1', 'site_2', 'site_3', 'site_4', 'site_5'.
setup = calibration_setup(file=ARPM$setup)
observed = calibration_data(setup=setup, path=ARPM$path)
#> Creating observed data list for calibration...
#> 
#> Loaded observed data for variables: 'site_1', 'site_2', 'site_3', 'site_4', 'site_5'.
forcing = as.matrix(read.csv(file.path(ARPM$path, "master", "environment.csv"), row.names=1))
control = list(maxit=20000, eps=sqrt(.Machine$double.eps), factr=sqrt(.Machine$double.eps))
```

Here, we also added a `control` list to increase the maximum number of
iterations for the algorithms and the tolerance for the convergence. Now
we can specify the `run_model()` function so we can simulate the model
from a parameter set and define the objective function using the
[`calibration_objFn()`](https://roliveros-ramos.github.io/calibrar/reference/calibration_objFn.md)
function:

``` r
run_model = function(par, forcing) {
  output = calibrar:::.PoissonMixedModel(par=par, forcing=forcing)
  output = c(output, list(gammas=par$gamma)) # adding gamma parameters for penalties
  return(output)
}
```

``` r
obj = calibration_objFn(model=run_model, setup=setup, observed=observed, forcing=forcing, aggregate=TRUE)
```

With these we can proceed to the parameter estimation. Here, we will
compare the performance of three BFGS type algorithms:

``` r
# real parameters
coef(ARPM)
#> $alpha
#> [1] 0.4
#> 
#> $beta
#> [1] -0.4
#> 
#> $gamma
#>  [1] -0.1252736456  0.0962670651  0.3390542168 -0.3522452588  0.0396026030
#>  [6]  0.0794698198  0.0058450990  0.5120546777  0.2514255423 -0.1069075371
#> [11] -0.1250454857  0.1827697374  0.2014399069  0.1438583647 -0.1209423322
#> [16]  0.1078108813 -0.0153661771  0.3699839121 -0.1709815103  0.0065274590
#> [21] -0.2050118962 -0.1964498151  0.0008203914 -0.0466854356 -0.0997776438
#> [26]  0.3099425926  0.0174993833  0.2637402261 -0.1962448238 -0.0491245175
#> [31] -0.2807867674  0.2881786295 -0.1962719982  0.2948489808 -0.1982394491
#> [36] -0.0188994689 -0.5750283368 -0.0493732202  0.0029488980 -0.3838175407
#> [41] -0.0575627487 -0.0693274891 -0.3679177169  0.1797177882 -0.2425710022
#> [46] -0.0437928463  0.1128536321 -0.1050868751  0.1488748450  0.0257963506
#> [51]  0.2976548514 -0.1325363902 -0.2321310003  0.0717548469 -0.0389692768
#> [56] -0.0590564016  0.0993280847  0.0969825576  0.0037569001  0.1269549129
#> [61]  0.1508888175  0.1667178070  0.1931522592  0.2587759936 -0.0273102042
#> [66] -0.0880277376 -0.2454567829 -0.0475306139 -0.1853716313  0.0822470712
#> [71] -0.0397729154 -0.1114899279 -0.1954314245  0.0121470497 -0.1194458959
#> [76] -0.2517897395 -0.2820094403  0.2202608159 -0.1354484053 -0.1524347015
#> [81] -0.0583497204 -0.1150377059 -0.0887883736 -0.0625140888 -0.1206008853
#> [86] -0.2187869440  0.1429412383 -0.0217624526 -0.2887595871  0.1612246527
#> [91] -0.3479670157 -0.0802640619 -0.0575164978 -0.1876814800  0.0575334278
#> [96] -0.3010802249  0.3038594027  0.0734818712  0.3399724818
#> 
#> $sd
#> [1] 0.2
#> 
#> $mu_ini
#> [1] 5.141664 4.158883 3.433987 4.430817 4.934474
```

``` r
lbfgsb1 = calibrate(par=ARPM$guess, fn=obj, method='L-BFGS-B', lower=ARPM$lower, upper=ARPM$upper, phases=ARPM$phase, control=control)
#> Using optimization method 'L-BFGS-B'.
#> Elapsed time: 1m 36.2s
#> Function value: -186095
#> Parameter values: 0.404 -0.415 -0.0827 0.128 0.254 -0.282 0.045 0.068 0.078 0.46 0.224 -0.0531 -0.102 0.148 0.301 0.0717 -0.0746 0.125 -0.0164 0.398 -0.188 0.0637 -0.213 -0.229 0.0825 -0.0407 -0.151 0.319 0.0263 0.294 -0.186 -0.0261 -0.297 0.276 -0.168 0.315 -0.216 -4.36e-05 -0.506 -0.0471 -0.0376 -0.276 -0.0826 -0.142 -0.304 0.0866 -0.162 -0.00317 0.131 -0.0579 0.0298 0.187 0.299 -0.0816 -0.347 0.112 -0.111 0.0628 0.107 0.0259 0.123 -0.0401 0.318 0.21 0.157 0.19 0.0758 -0.12 -0.0982 -0.136 -0.257 0.104 -0.0217 0.014 -0.207 -0.0266 -0.141 -0.159 -0.203 -0.00508 0.0131 -0.188 -0.048 -0.0728 -0.159 -0.112 -0.00272 -0.146 0.21 -0.177 -0.219 0.191 -0.259 -0.334 -0.0106 -0.174 0.0483 0.00983 0.235 -0.0243 0.414
#> 
#> Status: CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
lbfgsb2 = calibrate(par=ARPM$guess, fn=obj, method='Rvmmin', lower=ARPM$lower, upper=ARPM$upper, phases=ARPM$phase, control=control)
#> Using optimization method 'Rvmmin'.
#> Elapsed time: 34.67s
#> Function value: -186095
#> Parameter values: 0.397 -0.415 -0.0756 0.135 0.261 -0.276 0.0523 0.0756 0.0849 0.467 0.231 -0.0461 -0.0946 0.155 0.308 0.0792 -0.0679 0.131 -0.00869 0.405 -0.181 0.07 -0.205 -0.221 0.0894 -0.0333 -0.145 0.326 0.0332 0.301 -0.179 -0.0194 -0.291 0.284 -0.161 0.323 -0.209 0.00854 -0.5 -0.0402 -0.0303 -0.27 -0.0754 -0.135 -0.295 0.0925 -0.152 -0.00192 0.139 -0.0502 0.0356 0.194 0.307 -0.0754 -0.339 0.117 -0.102 0.0691 0.114 0.034 0.129 -0.0328 0.325 0.217 0.164 0.197 0.0826 -0.113 -0.0907 -0.128 -0.25 0.111 -0.0151 0.0209 -0.201 -0.0207 -0.134 -0.154 -0.194 0.000292 0.0197 -0.179 -0.0424 -0.0658 -0.151 -0.103 -0.00237 -0.129 0.209 -0.167 -0.206 0.187 -0.247 -0.321 -0.012 -0.164 0.0532 0.018 0.24 -0.00799 0.41
#> 
#> Status: Rvmminb appears to have converged
lbfgsb3 = calibrate(par=ARPM$guess, fn=obj, method='LBFGSB3', lower=ARPM$lower, upper=ARPM$upper, phases=ARPM$phase, control=control)
#> Using optimization method 'LBFGSB3'.
#>  Positive dir derivative in projection 
#>  Using the backtracking step
#> Elapsed time: 57.45s
#> Function value: -186093
#> Parameter values: 0.423 -0.415 -0.104 0.111 0.234 -0.303 0.0217 0.055 0.0675 0.44 0.203 -0.0732 -0.122 0.125 0.283 0.0572 -0.0945 0.109 -0.0361 0.382 -0.209 0.0423 -0.229 -0.249 0.0676 -0.0613 -0.174 0.309 -8.38e-05 0.28 -0.209 -0.0429 -0.326 0.264 -0.194 0.306 -0.241 -0.0119 -0.528 -0.0706 -0.0538 -0.296 -0.107 -0.158 -0.328 0.0917 -0.193 -0.0269 0.119 -0.0802 0.0183 0.168 0.282 -0.0984 -0.381 0.109 -0.144 0.046 0.0918 0.00202 0.112 -0.0693 0.307 0.192 0.145 0.169 0.0532 -0.136 -0.119 -0.154 -0.277 0.0881 -0.0462 0.00337 -0.231 -0.0444 -0.165 -0.173 -0.226 -0.0178 -0.00252 -0.223 -0.0665 -0.0912 -0.184 -0.13 -0.0129 -0.168 0.211 -0.214 -0.254 0.193 -0.282 -0.355 -0.02 -0.183 0.0314 -0.00983 0.208 -0.0465 0.414
#> 
#> Status: CONVERGENCE: Parameters differences below xtol
```

``` r
summary(ARPM, lbfgsb1, lbfgsb2, lbfgsb3, show_par = 1:3)
#>           method elapsed   value   fn   gr alpha   beta  gamma1
#> ARPM        data      NA -186178   NA   NA 0.400 -0.400 -0.1253
#> lbfgsb1 L-BFGS-B    96.2 -186095 2933 2933 0.404 -0.415 -0.0827
#> lbfgsb2   Rvmmin    34.7 -186095 4011  637 0.397 -0.415 -0.0756
#> lbfgsb3  LBFGSB3    57.4 -186093 1735 1735 0.423 -0.415 -0.1037
```

In this case, the best solution was found using the ‘Rvmmin’ algorithm,
which was also the faster (12.2s). Now, we can try to carry out the
parameter estimation in two phases:

``` r
phases = ARPM$phase
phases$gamma[] = 2
```

And re-do every optimization:

``` r
lbfgsb1p = calibrate(par=ARPM$guess, fn=obj, method='L-BFGS-B', lower=ARPM$lower, upper=ARPM$upper, phases=phases, control=control)
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 107 parameters are currently active.
#>  Using optimization method 'L-BFGS-B'.
#>  Phase 1 finished (0.08s)
#>  Function value: -169005
#>  Parameter values: 0.292 -0.295
#> 
#> - Phase 2: 101 out of 107 parameters are currently active.
#>  Using optimization method 'L-BFGS-B'.
#>  Phase 2 finished (1m 27.5s)
#>  Function value: -186095
#>  Parameter values: 0.4 -0.415 -0.078 0.133 0.258 -0.278 0.05 0.0733 0.0823 0.464 0.228 -0.0482 -0.0968 0.152 0.305 0.0767 -0.0702 0.129 -0.011 0.402 -0.183 0.0679 -0.208 -0.224 0.087 -0.0356 -0.147 0.323 0.0307 0.298 -0.181 -0.0216 -0.293 0.282 -0.163 0.32 -0.212 0.0064 -0.503 -0.0426 -0.0328 -0.272 -0.078 -0.137 -0.298 0.0917 -0.157 3.16e-05 0.136 -0.0533 0.0337 0.192 0.304 -0.0771 -0.342 0.117 -0.105 0.0674 0.112 0.0311 0.127 -0.0357 0.323 0.215 0.162 0.194 0.0803 -0.115 -0.0931 -0.131 -0.253 0.109 -0.0176 0.0188 -0.203 -0.0225 -0.137 -0.155 -0.198 -0.00116 0.0183 -0.182 -0.044 -0.0676 -0.154 -0.108 0.0013 -0.133 0.207 -0.17 -0.207 0.185 -0.25 -0.322 -0.0193 -0.16 0.0463 0.0206 0.228 -6.59e-05 0.403
#> 
#> Status: CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
lbfgsb2p = calibrate(par=ARPM$guess, fn=obj, method='Rvmmin', lower=ARPM$lower, upper=ARPM$upper, phases=phases, control=control)
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 107 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 1 finished (0.03s)
#>  Function value: -169005
#>  Parameter values: 0.292 -0.295
#> 
#> - Phase 2: 101 out of 107 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 2 finished (22.77s)
#>  Function value: -186095
#>  Parameter values: 0.396 -0.415 -0.0745 0.136 0.262 -0.275 0.0534 0.0766 0.086 0.468 0.232 -0.0449 -0.0937 0.156 0.309 0.0803 -0.0669 0.132 -0.00748 0.406 -0.18 0.071 -0.204 -0.22 0.0903 -0.0322 -0.144 0.327 0.0341 0.302 -0.178 -0.0183 -0.29 0.285 -0.16 0.324 -0.208 0.00981 -0.499 -0.0392 -0.0292 -0.269 -0.0744 -0.134 -0.294 0.0933 -0.153 0.0039 0.138 -0.0494 0.0369 0.195 0.308 -0.0742 -0.338 0.118 -0.1 0.0701 0.115 0.0354 0.13 -0.0315 0.325 0.218 0.165 0.198 0.0836 -0.112 -0.0896 -0.127 -0.249 0.112 -0.0137 0.0219 -0.2 -0.0195 -0.133 -0.152 -0.194 0.00193 0.021 -0.178 -0.0414 -0.0647 -0.15 -0.102 -0.000635 -0.128 0.21 -0.166 -0.205 0.188 -0.246 -0.321 -0.0105 -0.164 0.0537 0.0197 0.239 -0.00219 0.409
#> 
#> Status: Rvmminb appears to have converged
lbfgsb3p = calibrate(par=ARPM$guess, fn=obj, method='LBFGSB3', lower=ARPM$lower, upper=ARPM$upper, phases=phases, control=control)
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 107 parameters are currently active.
#>  Using optimization method 'LBFGSB3'.
#>  Phase 1 finished (0.10s)
#>  Function value: -169005
#>  Parameter values: 0.292 -0.295
#> 
#> - Phase 2: 101 out of 107 parameters are currently active.
#>  Using optimization method 'LBFGSB3'.
#>  Positive dir derivative in projection 
#>  Using the backtracking step 
#>  Positive dir derivative in projection 
#>  Using the backtracking step 
#>  Positive dir derivative in projection 
#>  Using the backtracking step
#>  Phase 2 finished (1m 3.0s)
#>  Function value: -186095
#>  Parameter values: 0.401 -0.415 -0.079 0.13 0.256 -0.279 0.0473 0.071 0.0825 0.463 0.227 -0.0499 -0.0995 0.151 0.304 0.0751 -0.0724 0.127 -0.0126 0.402 -0.185 0.0663 -0.209 -0.225 0.0858 -0.0375 -0.148 0.321 0.0295 0.297 -0.182 -0.0234 -0.295 0.28 -0.163 0.318 -0.215 0.00718 -0.505 -0.0468 -0.0326 -0.273 -0.0779 -0.14 -0.303 0.0963 -0.163 -0.00418 0.14 -0.0564 0.0306 0.191 0.306 -0.0797 -0.348 0.122 -0.11 0.0641 0.112 0.0302 0.128 -0.041 0.323 0.212 0.161 0.193 0.0788 -0.117 -0.0957 -0.131 -0.257 0.11 -0.0212 0.0199 -0.208 -0.0215 -0.14 -0.156 -0.202 -0.000204 0.0161 -0.185 -0.0462 -0.0672 -0.155 -0.11 -0.00264 -0.131 0.203 -0.177 -0.206 0.181 -0.244 -0.322 -0.0296 -0.156 0.0412 0.0192 0.22 0.0108 0.394
#> 
#> Status: CONVERGENCE: Parameters differences below xtol
```

``` r
summary(ARPM, lbfgsb1, lbfgsb1p, lbfgsb2, lbfgsb2p, lbfgsb3, lbfgsb3p, show_par = 1:3)
#>            method elapsed   value   fn   gr alpha   beta  gamma1
#> ARPM         data      NA -186178   NA   NA 0.400 -0.400 -0.1253
#> lbfgsb1  L-BFGS-B    96.2 -186095 2933 2933 0.404 -0.415 -0.0827
#> lbfgsb1p L-BFGS-B    87.6 -186095 2695 2695 0.400 -0.415 -0.0780
#> lbfgsb2    Rvmmin    34.7 -186095 4011  637 0.397 -0.415 -0.0756
#> lbfgsb2p   Rvmmin    22.8 -186095 2650  427 0.396 -0.415 -0.0745
#> lbfgsb3   LBFGSB3    57.4 -186093 1735 1735 0.423 -0.415 -0.1037
#> lbfgsb3p  LBFGSB3    63.1 -186095 1964 1964 0.401 -0.415 -0.0790
```

For all the algorithms, we can see and improvement in the solution found
and even a speed up in the time needed for the optimization for the
first two methods (‘L-BFGS-B’, ‘Rvmmin’). Finally, we can try to carry
out the optimization in two phases, but using different algorithms for
each phase. In principle, some algorithms may be faster or find a better
initial starting point for the final search.

``` r
mix1 = calibrate(par=ARPM$guess, fn=obj, method=c('hjn', 'Rvmmin'), lower=ARPM$lower, upper=ARPM$upper, phases=phases, control=control)
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 107 parameters are currently active.
#>  Using optimization method 'hjn'.
#>  Phase 1 finished (1.26s)
#>  Function value: -169005
#>  Parameter values: 0.294 -0.298
#> 
#> - Phase 2: 101 out of 107 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 2 finished (33.22s)
#>  Function value: -186095
#>  Parameter values: 0.397 -0.415 -0.0756 0.135 0.261 -0.276 0.0523 0.0757 0.085 0.467 0.231 -0.0459 -0.0946 0.155 0.308 0.0793 -0.0678 0.131 -0.0087 0.405 -0.181 0.0701 -0.205 -0.221 0.0894 -0.0333 -0.145 0.326 0.0331 0.301 -0.179 -0.0194 -0.291 0.284 -0.161 0.323 -0.209 0.00865 -0.5 -0.0401 -0.0303 -0.27 -0.0754 -0.135 -0.295 0.0923 -0.153 0.000349 0.139 -0.0503 0.0357 0.194 0.307 -0.0753 -0.339 0.117 -0.101 0.069 0.114 0.0343 0.129 -0.033 0.325 0.217 0.164 0.197 0.0826 -0.113 -0.0907 -0.128 -0.25 0.111 -0.015 0.0208 -0.201 -0.0205 -0.134 -0.154 -0.195 0.00193 0.0189 -0.178 -0.0424 -0.0657 -0.151 -0.103 -0.00121 -0.129 0.209 -0.167 -0.205 0.187 -0.247 -0.321 -0.0125 -0.164 0.0532 0.0177 0.24 -0.00831 0.411
#> 
#> Status: Rvmminb appears to have converged
mix2 = calibrate(par=ARPM$guess, fn=obj, method=c('Nelder-Mead', 'Rvmmin'), lower=ARPM$lower, upper=ARPM$upper, phases=phases, control=control)
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 107 parameters are currently active.
#>  Using optimization method 'Nelder-Mead'.
#>  Phase 1 finished (0.01s)
#>  Function value: -124445
#>  Parameter values: 0.308 -0.358
#> 
#> - Phase 2: 101 out of 107 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 2 finished (7.73s)
#>  Function value: -186095
#>  Parameter values: 0.399 -0.415 -0.0775 0.133 0.258 -0.277 0.0504 0.0735 0.083 0.465 0.229 -0.0479 -0.0964 0.153 0.306 0.0775 -0.0697 0.129 -0.0104 0.403 -0.182 0.0681 -0.207 -0.223 0.0874 -0.0351 -0.146 0.324 0.03 0.3 -0.18 -0.0213 -0.293 0.282 -0.163 0.321 -0.211 0.00609 -0.502 -0.0418 -0.0322 -0.271 -0.0771 -0.137 -0.296 0.0894 -0.156 0.00341 0.135 -0.0515 0.0305 0.194 0.306 -0.0774 -0.341 0.116 -0.103 0.0671 0.114 0.0282 0.129 -0.0344 0.323 0.215 0.162 0.195 0.0804 -0.114 -0.0933 -0.13 -0.251 0.109 -0.0147 0.0147 -0.2 -0.023 -0.136 -0.155 -0.194 0.00448 -2.57e-05 -0.173 -0.0452 -0.0652 -0.153 -0.105 -0.00109 -0.135 0.208 -0.167 -0.208 0.186 -0.249 -0.324 -0.0105 -0.166 0.0554 0.00559 0.24 -0.00624 0.41
#> 
#> Status: Rvmminb appears to have converged
mix3 = calibrate(par=ARPM$guess, fn=obj, method=c('CG', 'Rvmmin'), lower=ARPM$lower, upper=ARPM$upper, phases=phases, control=control)
#> Parameter estimation in two phases.
#> 
#> - Phase 1: 2 out of 107 parameters are currently active.
#>  Using optimization method 'CG'.
#>  Phase 1 finished (0.03s)
#>  Function value: -1824.73
#>  Parameter values: 2.33e+03 2.36e+03
#> 
#> - Phase 2: 101 out of 107 parameters are currently active.
#>  Using optimization method 'Rvmmin'.
#>  Phase 2 finished (55.75s)
#>  Function value: -186095
#>  Parameter values: 0.397 -0.415 -0.0756 0.135 0.261 -0.276 0.0523 0.0756 0.085 0.467 0.231 -0.0461 -0.0946 0.155 0.308 0.0793 -0.0679 0.131 -0.00868 0.405 -0.181 0.07 -0.205 -0.221 0.0894 -0.0333 -0.145 0.326 0.0331 0.301 -0.179 -0.0195 -0.291 0.284 -0.161 0.323 -0.209 0.00867 -0.5 -0.0402 -0.0304 -0.27 -0.0755 -0.135 -0.295 0.0924 -0.153 0.000323 0.138 -0.0502 0.0355 0.194 0.307 -0.0753 -0.339 0.117 -0.101 0.069 0.114 0.034 0.129 -0.0329 0.325 0.217 0.164 0.197 0.0826 -0.113 -0.0907 -0.128 -0.25 0.111 -0.0153 0.021 -0.201 -0.0206 -0.134 -0.153 -0.194 0.000435 0.0196 -0.179 -0.0424 -0.0658 -0.151 -0.103 -0.0018 -0.129 0.209 -0.167 -0.206 0.187 -0.247 -0.321 -0.013 -0.164 0.0534 0.0176 0.24 -0.00723 0.41
#> 
#> Status: Rvmminb appears to have converged
```

``` r
summary(ARPM, lbfgsb2, lbfgsb2p, mix1, mix2, mix3, show_par = 1:3)
#>          method elapsed   value   fn   gr alpha   beta  gamma1
#> ARPM       data      NA -186178   NA   NA 0.400 -0.400 -0.1253
#> lbfgsb2  Rvmmin   34.67 -186095 4011  637 0.397 -0.415 -0.0756
#> lbfgsb2p Rvmmin   22.81 -186095 2650  427 0.396 -0.415 -0.0745
#> mix1     Rvmmin   34.49 -186095 2498  635 0.397 -0.415 -0.0756
#> mix2     Rvmmin    7.74 -186095  410  149 0.399 -0.415 -0.0775
#> mix3     Rvmmin   55.78 -186095 4011 1070 0.397 -0.415 -0.0756
```

Here, we can see that every combination find essentially the same
solution, but the combination using the ‘CG’ (conjugated gradient) first
required less function and gradient evaluations, being faster.

Please, refer to `vignette(package="calibrar")` for additional vignettes
or to the [calibrar
website](https://roliveros-ramos.github.io/calibrar/) for more details.
