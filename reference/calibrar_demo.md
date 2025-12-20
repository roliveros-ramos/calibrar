# Demos for the calibrar package

Creates demo files able to be processed for a full calibration using the
calibrar package

## Usage

``` r
calibrar_demo(path = NULL, model = NULL, ...)
```

## Arguments

- path:

  Path to create the demo files

- model:

  Model to be used in the demo files, see details.

- ...:

  Additional parameters to be used in the construction of the demo
  files.

## Value

A list with the following elements:

- path:

  Path were the files were saved

- par:

  Real value of the parameters used in the demo

- setup:

  Path to the calibration setup file

- guess:

  Values to be provided as initial guess to the calibrate function

- lower:

  Values to be provided as lower bounds to the calibrate function

- upper:

  Values to be provided as upper bounds to the calibrate function

- phase:

  Values to be provided as phases to the calibrate function

- constants:

  Constants used in the demo, any other variable not listed here.

- value:

  NA, set for compatibility with summary methods.

- time:

  NA, set for compatibility with summary methods.

- counts:

  NA, set for compatibility with summary methods.

## Details

Current implemented models are:

- PoissonMixedModel:

  Poisson Autoregressive Mixed model for the dynamics of a population in
  different sites: \$\$log(\mu\_{i, t+1}) = log(\mu\_{i, t}) + \alpha +
  \beta X\_{i, t} + \gamma_t\$\$ where \\\mu\_{i, t}\\ is the size of
  the population in site \\i\\ at year \\t\\, \\X\_{i, t}\\ is the value
  of an environmental variable in site \\i\\ at year \\t\\. The
  parameters to estimate were \\\alpha\\, \\\beta\\, and \\\gamma_t\\,
  the random effects for each year, \\\gamma_t \sim N(0,\sigma^2)\\, and
  the initial population at each site \\\mu\_{i, 0}\\. We assumed that
  the observations \\N\_{i,t}\\ follow a Poisson distribution with mean
  \\\mu\_{i, t}\\.

- PredatorPrey:

  Lotka Volterra Predator-Prey model. The model is defined by a system
  of ordinary differential equations for the abundance of prey \$N\$ and
  predator \$P\$: \$\$\frac{dN}{dt} = rN(1-N/K)-\alpha NP\$\$
  \$\$\frac{dP}{dt} = -lP + \gamma\alpha NP\$\$ The parameters to
  estimate are the prey’s growth rate \\r\\, the predator’s mortality
  rate \\l\\, the carrying capacity of the prey \\K\\ and \\\alpha\\ and
  \\\gamma\\ for the predation interaction. Uses `deSolve` package for
  numerical solution of the ODE system.

- SIR:

  Susceptible-Infected-Recovered epidemiological model. The model is
  defined by a system of ordinary differential equations for the number
  of susceptible \$S\$, infected \$I\$ and recovered \$R\$ individuals:
  \$\$\frac{dS}{dt} = -\beta S I/N\$\$ \$\$\frac{dI}{dt} = \beta S I/N
  -\gamma I\$\$ \$\$\frac{dR}{dt} = \gamma I\$\$ The parameters to
  estimate are the average number of contacts per person per time
  \\\beta\\ and the instant probability of an infectious individual
  recovering \\\gamma\\. Uses `deSolve` package for numerical solution
  of the ODE system.

- IBMLotkaVolterra:

  Stochastic Individual Based Model for Lotka-Volterra model. Uses `ibm`
  package for the simulation.

## References

Oliveros-Ramos and Shin (2014)

## Author

Ricardo Oliveros–Ramos

## Examples

``` r
if (FALSE) { # \dontrun{

summary(ahr)
set.seed(880820)
path = NULL # NULL to use the current directory
# create the demonstration files
demo = calibrar_demo(path=path, model="PredatorPrey", T=100) 
# get calibration information
calibration_settings = calibration_setup(file = demo$setup)
# get observed data
observed = calibration_data(setup = calibration_settings, path=demo$path)
# Defining 'run_model' function
run_model = calibrar:::.PredatorPreyModel
# real parameters
cat("Real parameters used to simulate data\n")
print(unlist(demo$par)) # parameters are in a list
# objective functions
obj  = calibration_objFn(model=run_model, setup=calibration_settings, observed=observed, T=demo$T)
obj2 = calibration_objFn(model=run_model, setup=calibration_settings, observed=observed, 
T=demo$T, aggregate=TRUE)
cat("Starting calibration...\n")
cat("Running optimization algorithms\n", "\t")
cat("Running optim AHR-ES\n")
ahr = calibrate(par=demo$guess, fn=obj, lower=demo$lower, upper=demo$upper, phases=demo$phase)
summary(ahr)
} # } 
```
