#!/usr/bin/env Rscript

# USER CONFIGURATION - MODIFY THIS SECTION --------------------------------

# Simple objective function (replace with your actual model)
objFn = function(x) {
  sum(x^2) + rnorm(1, sd=0.01)  # add small noise for replicates
}

# Parameter setup
par_guess = rep(0.5, 3)
par_min = rep(-5, 3)
par_max = rep(5, 3)
par_phase = c(1, 1, 2)

method = NULL # use default method

control = list() # additional arguments

# DO NOT MODIFY BELOW THIS LINE -------------------------------------------

library(calibrar)
library(parallel)

if(!exists(".args", mode = "character")) .args = commandArgs(trailingOnly=TRUE)
# Get command line arguments
ncores = .get_command_argument(.args, "ncores", default=1)
replicates = .get_command_argument(.args, "replicates", default=1)
if(replicates > 1) {
  message("Starting calibration with ", ncores, " cores, ", replicates, " replicates.")
} else {
  message("Starting calibration with ", ncores, " cores.")
}

if(!exists("control")) control = list()
control$ncores = ncores

# Setup parallel backend
if(ncores > 1) {
  cl = makeCluster(ncores)
  # Export necessary objects to cluster
  clusterExport(cl, "objFn")
  # Load calibrar on workers
  clusterEvalQ(cl, library(calibrar))
} else {
  cl = NULL
}

# Run calibration
result = calibrate(
  par = par_guess,
  fn = objFn,
  method = method,
  lower = par_min,
  upper = par_max,
  phases = par_phase,
  replicates = replicates,
  control = control,
  parallel = (ncores > 1)
)

# cleanup
if(!is.null(cl)) stopCluster(cl)

message("Calibration completed.")
print(result)
