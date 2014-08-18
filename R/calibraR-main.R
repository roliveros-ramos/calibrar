# TO_DO: replicates for fn
# multiple phases
# check par: length(par)==1 provide the number of parameters, guess=rep(NA, par)?
# add npar? Check compatibility or use it to set the calibration, guess=rep(NA, npar)


optimES = function (par, fn, gr = NULL, ..., method = "default", 
                    lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, restart=NULL) {
  
  npar = length(par)

  active = .checkActive(active=active, npar=npar)
  bounds = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess  = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  
  isActive = which(active)
  
  par    = guess[isActive]
  lower  = bounds$lower[isActive]
  upper  = bounds$upper[isActive]
  
  npar = length(par)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    fn(parx, ...)
  }
  
  control = .checkControl(control=control, method=method, par=par, fn=fn1, active=active)
  
  if(control$REPORT>0) {
    trace = list()
    trace$control = control
    trace$value = numeric(control$maxgen)
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
    if(control$REPORT>0) trace$value[opt$gen] = control$aggFn(fn1(opt$MU), control$weights)
    
  }
  
  value = control$aggFn(fn1(opt$MU), control$weights)
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)
  
  output = list(par=opt$MU, value=value, counts=opt$counts, 
                trace=trace)
  
  return(output)
  
}


calibrate = function(par, fn, ..., aggFn = NULL, method = "default",
                     replicates=1, lower = -Inf, upper = Inf, phases=NULL, 
                     gr = NULL, control = list(), hessian = FALSE, 
                     restart=NULL) {
  
  npar = length(par)
  
  fn = match.fun(fn)
  
  phases     = .checkPhases(phases=phases, npar=npar)
  bounds     = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess      = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)

  par    = guess
  lower  = bounds$lower
  upper  = bounds$upper
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
    
  }
  
  output = c(temp, output)
  class(output) = c("calibrar", class(output))
  
  return(output)
  
}






