
optimEA = function (par, fn, gr = NULL, ..., method = "default", 
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
    
  }
  
  value = control$aggFn(fn1(opt$MU), control$weights)
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)
  
  output = list(par=opt$MU, value=value, counts=opt$counts)
  
  return(output)
  
}



