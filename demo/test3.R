require(calibrar)
require(optimx)
require(cmaes)

set.seed(12345) # updated to T=20 and L=40 for comparative purposes.

path = "private" # path where to create the demo files

# create the demostration files
# demo = calibrarDemo(path=path, model="PoissonMixedModel", L=4, T=20) 

# set.seed(12345) # updated to T=20 and L=40 for comparative purposes.
# Parameter's information
parInfo = read.csv(file.path(demo$path, "parInfo.csv"), row.names=1)

# get calibration information
calibrationInfo = getCalibrationInfo(path=demo$path)

# get observed data
observed = getObservedData(info=calibrationInfo, path=demo$path)

# read forcings for the model
forcing = read.csv(file.path(demo$path, "master", "environment.csv"),
                   row.names=1)

# Defining 'runModel' function
  runModel = function(par, forcing) {
  
    # forcing is a matrix with the values of environmental variables
    T = nrow(forcing) # get number of time steps
    L = ncol(forcing) # get number of sites
    
    # create parameter list in a format readable by model
    parList = list()
    parList$alpha = par[1]
    parList$beta  = par[2]
    parList$gamma = par[2 + seq_len(T-1)]
    parList$sd = par[T+2]
    parList$mu_ini = par[T + 2 + seq_len(L)]
  
    # get the model
    model = calibrar:::.PoissonMixedModel
    # run the model
    output = model(par=parList, forcing=forcing)
    output = as.list(as.data.frame(output)) # return a list with the results
    # names of the outputs matching observed data names
    names(output) = paste0("site_", sprintf(paste0("%0", ceiling(log10(L+1)), "d"), seq_len(L)))
    
    output = c(output, list(gammas=parList$gamma))
    return(output)
    
  }

x = runModel(parInfo$guess, forcing) # results from first guess
# the runModel function produces a list with the simulated outputs
# matching the names of the observed variables
# par(mfrow=c(2,2))
print(x)
# lapply(x, plot, type="l")

names(x)

obj = createObjectiveFunction(runModel=runModel, info=calibrationInfo, observed=observed,
                              forcing=forcing)
obj2 = createObjectiveFunction(runModel=runModel, info=calibrationInfo, observed=observed,
                              forcing=forcing, aggregate=TRUE)
cat("Starting calibration...\n")

calib = calibrate(par=parInfo$guess, fn=obj, lower=parInfo$lower, 
                  upper=parInfo$upper, phases=parInfo$phase,
                  control=list(weights=calibrationInfo$weights, REPORT=1, trace=1))

# calib2 = optimES(par=parInfo$guess, fn=obj, lower=parInfo$lower, 
#                   upper=parInfo$upper, control=list(weights=calibrationInfo$weights))
# 
# calib3 = optimx(par=parInfo$guess, fn=obj2, lower=parInfo$lower, 
#                  upper=parInfo$upper, method="L-BFGS-B")
# 
# calib4 = cma_es(par=parInfo$guess, fn=obj2, lower=parInfo$lower, 
#                 upper=parInfo$upper)
# 
# calib5 = optim(par=parInfo$guess, fn=obj2)
# 
# calib6 = optim(par=parInfo$guess, fn=obj2, method="SANN")

print(calib)




