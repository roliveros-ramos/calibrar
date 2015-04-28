.calibrar = function (par, fn, gr = NULL, ..., lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, method = "default") {
  
  
  npar = length(par)
 
  # check active parameters
  active = .checkActive(active=active, npar=npar)
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  # update to active parameters only
  guess  = par
  par    = guess[isActive]
  lower  = lower[isActive]
  upper  = upper[isActive]
  
  npar = length(par)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  control = .checkControl(control=control, method=method, par=par, fn=fn, active=active)
  
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    fn(parx, ...)/control$fnscale
  }

  output = .optimES(par=par, fn=fn1, lower=lower, upper=upper, control=control)

  # reshaping full parameters
  paropt = guess
  paropt[isActive] = output$ppar 
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  # final outputs
  output = list(par=paropt, output, active=list(par=isActive, flag=activeFlag))
  
  class(output) = c("optimES.result", class(output)) # change
  
  return(output)
  
}


# optimES internal --------------------------------------------------------

.optimES = function(par, fn, lower, upper, control, hessian=FALSE) {
  
  # get restart for the current phase
  restart = .restartCalibration(control) # flag: TRUE or FALSE
  if(isTRUE(restart)) {
    
    res = .getRestart(control=control)
    opt   = res$opt
    trace = res$trace
    
  } else {
    
    opt = .newOpt(par=par, lower=lower, upper=upper, control=control)
    
    trace = NULL
    
    if(control$REPORT>0 & control$trace>0) {
      
      trace = list()
      trace$control = control
      trace$par = matrix(NA, nrow=control$maxgen, ncol=length(isActive))
      trace$value = rep(NA, control$maxgen)
      trace$best  = rep(NA, control$maxgen)
      
      if(control$trace>1) {
        trace$sd = matrix(NA, nrow=control$maxgen, ncol=length(isActive))   
        trace$step = rep(NA, control$maxgen)     
      }
      
      if(control$trace>2) trace$opt = vector("list", control$maxgen)
      
    } 
    
  } 
  
  # start new optimization
  while(isTRUE(.continueEvolution(opt, control))) {
    
    opt$gen  = opt$gen + 1
    opt$ages = opt$ages + 1
    
    # create a new population
    
    if(all(opt$SIGMA==0)) break
    opt$pop = .createPopulation(opt)
    
    # evaluate the function in the population: evaluate fn, aggregate fitness
    
    opt$fitness = .calculateFitness(opt, fn=fn)
    
    # select best 'individuals'
    
    opt$selected = .selection(opt)
    
    # create the new parents: MU and SD
    
    opt = .calculateOptimalWeights(opt)
    opt = .updatePopulation(opt)
    
    # save detailed outputs
    if(control$REPORT>0 & control$trace>0) {
      
      trace$par[opt$gen, ] = opt$MU
      trace$best[opt$gen]  = opt$selected$best$fit.global
      
      if(control$trace>1) {
        trace$sd[opt$gen, ]  = opt$SIGMA
        trace$step[opt$gen]  = opt$step       
      }
      
      if(opt$gen%%control$REPORT==0) {
        trace$value[opt$gen] = control$aggFn(fn1(opt$MU), control$weights)
        if(control$trace>2) trace$opt[[opt$gen]] = opt
      }
      
    }
    
    # save restart
    .createRestartFile(opt=opt, trace=trace, control=control)
    
    if(control$verbose & opt$gen%%control$REPORT==0) 
      .messageByGen(opt, trace)
    
  } # end generations loop
  
  partial = fn(opt$MU)
  value = control$aggFn(x=partial, w=control$weights) # check if necessary
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)
  
  output = list(ppar=opt$MU, value=value, counts=opt$counts, 
                trace=trace, partial=partial)
  
  return(output)
  
}
