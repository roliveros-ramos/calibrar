
optimEA = function (par, fn, gr = NULL, ..., method = "default", 
                    lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, restart=NULL) {
  
  fn   = match.fun(fn)
  npar = length(par)
  
  bounds = .checkBounds(lower=lower, upper=upper, npar=npar)
  
  lower  = bounds$lower
  upper  = bounds$upper
  par    = .checkOpt(par=par, lower=lower, upper=upper)
  
  control = .checkControl(control=control, method=method, par=par, fn=fn, ...)
  
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
    
    opt$fitness = .calculateFitness(opt, fn=fn, ...)
    
    # select best 'individuals'
    
    opt$selected = .selection(opt)
        
    # create the new parents: MU and SD
    
    opt = .calculateOptimalWeights(opt)
    opt = .updatePopulation(opt)
    
    # save status of the population (for restart)
    
    # save detailed outputs
    
  }
  
  value = control$aggFn(fn(opt$MU), control$weights)
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)
  
  output = list(par=opt$MU, value=value, counts=opt$counts)
  
  return(output)
  
}



