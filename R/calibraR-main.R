# TO_DO: replicates for fn
# check par: length(par)==1 provide the number of parameters, guess=rep(NA, par)?
# add npar? Check compatibility or use it to set the calibration, guess=rep(NA, npar)
# Restart
# Should function be taken from restart or arguments?
# Should we take control from new control or restart?
# Now: to reproduce actual functioning, think on mutating objects.

# save restart in run dir if not null.
# use default restart.file if not defined, build over model? parse function name?
# this has to be inherited from calibrar() to include phase and model
# define saveRestart?
# messages by REPORT

# model-specific packages can be designed over calibrar.



optimES = function (par, fn, gr = NULL, ..., method = "default", 
                    lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, restart=NULL) {
  
  npar = length(par)

  if(is.null(names(par))) names(par) = .printSeq(npar, preffix="par")
  
  active = .checkActive(active=active, npar=npar)
  bounds = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess  = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  par    = guess[isActive]
  lower  = bounds$lower[isActive]
  upper  = bounds$upper[isActive]
  
  npar = length(par)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  
  control = .checkControl(control=control, method=method, par=par, fn=fn, active=active)
  
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    fn(parx, ...)/control$fnscale
  }
  
  
  if(control$REPORT>0 & control$trace>0) {
    trace = list()
    trace$control = control
    trace$value = rep(NA, control$maxgen)
    trace$step = rep(NA, control$maxgen)
    trace$par = matrix(NA, nrow=control$maxgen, ncol=length(isActive))
    trace$sd = matrix(NA, nrow=control$maxgen, ncol=length(isActive))
    trace$opt = vector("list", control$maxgen)
  } else trace=NULL

  
  # opt = get restart, a method for a file (character) or a restart class
  
  opt = if(!is.null(restart)) # TO_DO
    .getRestart(restart=restart) else 
      .newOpt(par=par, lower=lower, upper=upper, control=control) # here par=NA?
  
  while(isTRUE(.continueEvolution(opt, control))) {
    
    opt$gen  = opt$gen + 1
    opt$ages = opt$ages + 1
    
    # create a new population
    
    if(all(opt$SIGMA==0)) break
    opt$pop = .createPopulation(opt)
    
    # evaluate the function in the population: evaluate fn, aggregate fitness
    
    opt$fitness = .calculateFitness(opt, fn=fn1)
    
    # select best 'individuals'
    
    opt$selected = .selection(opt)
        
    # create the new parents: MU and SD
    
    opt = .calculateOptimalWeights(opt)
    opt = .updatePopulation(opt)
    
    # save status of the population (for restart)
    
    # save detailed outputs
    if(control$REPORT>0 & control$trace>0) {
      
      if(opt$gen%%control$REPORT==0) {
        trace$value[opt$gen] = control$aggFn(fn1(opt$MU), control$weights)
        trace$step[opt$gen] = opt$step
        trace$par[opt$gen, ] = opt$MU
        trace$sd[opt$gen, ] = opt$SIGMA
        trace$opt[[opt$gen]] = opt
        .messageByGen(opt, trace)
      }
      
    }
    
  }
  
  value = control$aggFn(x=fn1(opt$MU), w=control$weights)
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)

  output = list(par=opt$MU, value=value, counts=opt$counts, 
                trace=trace, partial=fn1(opt$MU), 
                active=list(par=isActive, flag=activeFlag))
  
  class(output) = c("optimES.result", class(output))
  
  return(output)
  
}


calibrate = function(par, fn, ..., aggFn = NULL, method = "default",
                     phases = NULL, replicates=1, lower = -Inf, upper = Inf,  
                     gr = NULL, control = list(), hessian = FALSE, 
                     restart = NULL) {
  
  npar = length(par)
  
  fn = match.fun(fn)
  
  phases     = .checkPhases(phases=phases, npar=npar)
  bounds     = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess      = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)

  par     = guess
  lower   = bounds$lower
  upper   = bounds$upper
  nphases = max(phases, na.rm=TRUE)

  replicates = .checkReplicates(replicates, nphases) 
  
  output = list()
  
  for(phase in seq_len(nphases)) {
      
    # call optimEA
    active = (phases <= phase) # NAs are corrected in optimES 
    temp = optimES(par=par, fn=fn, gr = NULL, ..., method = method, 
                   lower = lower, upper = upper, active=active, 
                   control = control, hessian = hessian, restart=restart)
   
    output$phases[[phase]] = temp # trim?
    par[which(active)] = temp$par
    control = .updateControl(control=control, opt=temp, method=method)  # update CVs? 

    cat(sprintf("Phase %d finished (%d of %d parameters active).\n",
                phase, sum(active, na.rm=TRUE), npar))
    print(temp$par)
  }
  
  isActive = !is.na(phases) & (phases>=1)
  paropt = guess
  paropt[isActive] = temp$par
  
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  newNames = rep("*", npar)
  newNames[isActive] = ""
  
  
  
  final = list(par=paropt, value=temp$value, counts=temp$counts, 
               trace=NULL, partial=temp$partial,
               active=isActive)
  
  output = c(final, output)
  class(output) = c("calibrar.result", class(output))
  
  return(output)
  
}


# getObservedData


getObservedData = function(info, path, data.folder="data", ...) {
  
  observed  = list()
  variables = info$variable
  
  useData       = as.logical(info$useData)
  
  cat("Creating observed data list for calibration...","\n")
  
  for(var in 1:nrow(info)) {
    
    cat(paste0("Variable: ", variables[var], "\n"))
    var.path        = file.path(path, data.folder, paste0(variables[var],".csv"))
    datos           = if(useData[var]) .read.csv3(var.path, ...) else NA
    observed[[var]] = datos
    
  }
  
  names(observed) = variables
  
  return(observed)
  
}

# get Calibration Info

getCalibrationInfo = function(path, file="calibrationInfo.csv", 
                              stringsAsFactors=FALSE, ...) {
  
  caliPath = file.path(path, file)
  calibrationInfo = read.csv(caliPath, stringsAsFactors=FALSE, ...)
  
  fullNames = c("variable", "type", "calibrate", "weights", "useData")  
  doesNotMatch = !(names(calibrationInfo) %in% fullNames)
  dnm = names(calibrationInfo)[doesNotMatch]
  
  isMissing = !(fullNames %in% names(calibrationInfo))
  im = fullNames[isMissing]
  
  sdnm = if(length(dnm)>1) " columns do " else " column does "
  sim  = if(length(im)>1) " variables are " else " variable is "
  msg1 = paste0("Error in ", caliPath, " file (", paste(sapply(dnm, sQuote), collapse=", "), 
                sdnm, "not match).")
  msg2 = paste0("Error in ", caliPath, " file (", paste(sapply(im, sQuote), collapse=", "), 
                sim, "missing).")
  
  if(any(doesNotMatch)) stop(msg1)
  if(any(isMissing)) stop(msg2)
  
  # cating correct data types
  calibrationInfo$variable  = as.character(calibrationInfo$variable)
  calibrationInfo$type      = as.character(calibrationInfo$type)
  calibrationInfo$calibrate = as.logical(calibrationInfo$calibrate)
  calibrationInfo$weights   = as.numeric(calibrationInfo$weights)
  calibrationInfo$useData   = as.logical(calibrationInfo$useData)
  
  return(calibrationInfo)
}

# createObjectiveFunction

createObjectiveFunction = function(runModel, info, observed, aggFn=.weighted.sum, 
                                   aggregate=FALSE, ...) {

  fn   = match.fun(runModel)
  aggFn = match.fun(aggFn)
  
  force(observed)
  force(info)
  force(aggregate)
  
  # check for names in observed and simulated
  
  fn1  = function(par) {
    aggFn = match.fun(aggFn)
    simulated = fn(par, ...)
    # apply fitness to all outputs
    output = .calculateObjetiveValue(obs=observed, sim=simulated, info=info)
    if(isTRUE(aggregate)) output = aggFn(x=output, w=info$weights)
    return(output)
  }
  
  return(fn1) 
  
}

.calculateObjetiveValue = function(obs, sim, info) {
  fit = NULL
  for(j in seq_len(nrow(info))) {
    if(!info$calibrate[j]) next
    var = info$variable[j]
    fit = c(fit, .fitness(obs=obs[[var]], sim=sim[[var]], FUN=info$type[j]))
  }
  names(fit) = info$variable[which(info$calibrate)]
  return(fit)
}

.fitness = function(obs, sim, FUN, ...) {
  FUN = match.fun(FUN)
  output = FUN(obs=obs, sim=sim, ...)
  return(output)
}

# calibrar Demo

calibrarDemo = function(path=NULL, model=NULL,  ...) {
  
  if(is.null(path)) path = getwd()
  if(is.null(model)) {
    model = "default"
    warning("Using default demo 'PoissonMixedModel'")
  }
  
  output = switch(model, 
                  Poisson = .generatePoissonMixedModel(path=path, seed=seed, ...),
                  .generatePoissonMixedModel(path=path, seed=seed, ...)  
  )
  return(output)                
  
}

