
# wrapper for cma_es ------------------------------------------------------

.cmaes = function(par, fn, lower, upper, control) {
  
  npar = length(unlist(par))
  
  if(!is.null(control$trace)) control$trace = ifelse(control$trace>1, TRUE, FALSE)
  
  # all default values!
  con = list(trace = FALSE, fnscale = 1, stopfitness = -Inf, maxit = 100 * npar^2,
             sigma = 0.5, keep.best = TRUE, vectorized = FALSE, diag = FALSE,
             lambda = 4 + floor(3 * log(npar)))
  
  con = within(con, {
    stop.tolx = 1e-12 * sigma
    diag.sigma = diag
    diag.eigen = diag
    diag.value = diag
    diag.pop = diag
    mu = floor(lambda/2)
    weights = log(mu + 1) - log(1:mu)
    weights = weights/sum(weights)
    mueff = sum(weights)^2/sum(weights^2)
    ccum = 4/(npar + 4)
    cs = (mueff + 2)/(npar + mueff + 3)
    ccov.mu = mueff
    ccov.1 = (1/ccov.mu) * 2/(npar + 1.4)^2 + (1 - 1/ccov.mu) * ((2 * ccov.mu - 1)/((npar + 2)^2 + 2 * ccov.mu))
    damps = 1 + 2 * max(0, sqrt((mueff - 1)/(npar + 1)) - 1) + cs
  })
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(cmaes::cma_es(par=par, fn=fn, lower=lower, upper=upper, control=control))
  
  if(is.null(output$par)) {
    output$par = relist(rep(NA, npar), skeleton=par)
    if(!is.finite(output$value)) warning("Infinite value reached for fn.")
  }
  
  return(output)
  
}



# wrapper for genSA -------------------------------------------------------

.genSA = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  con = list(maxit = 5000, threshold.stop = NULL, temperature = 5230, 
             visiting.param = 2.62, acceptance.param = -5, max.time = NULL, 
             nb.stop.improvement = 1e+06, smooth = TRUE, max.call = 1e+07, 
             simple.function = FALSE, trace.fn = NULL, verbose = FALSE, 
             trace.mat = TRUE, seed = -100377, high.dim = TRUE, tem.restart = 0.1,
             markov.length = 2 * length(lower))
  
  control = check_control(control=control, default=con)
  
  xoutput = suppressWarnings(GenSA::GenSA(par=par, fn=fn, lower=lower, 
                                         upper=upper, control=control))
  
  output = list()
  output$par  = xoutput$par
  output$value = xoutput$value
  output$counts = c('function'=xoutput$counts, 'gradient'=0)	
  output$convergence = 0
  output$message = NA
  output$hessian = NULL
  
  return(output)
  
}



# wrapper for DEoptim -----------------------------------------------------

.DE = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(isTRUE(control$parallel)) control$parallelType="auto"
  
  con = DEoptim::DEoptim.control()
  control = check_control(control=control, default=con)
  
  # we pass 'par' as initial population (gaussian with uniform-like variance).
  if(is.null(control$initialpop)) {
    control$initialpop = t(.createRandomPopulation(n=control$NP, mean=par, lower=lower, upper=upper)) 
  }

  # think how to pass 'par'. Initial pop?
  xoutput = suppressWarnings(DEoptim::DEoptim(fn=fn, lower=lower, upper=upper, control=control))
  
  output = list()
  output$par  = xoutput$optim$bestmem
  output$value = xoutput$optim$bestval
  output$counts = c('function'=xoutput$optim$nfeval,
                    gradient=NA)
  output$convergence = NA
  output$message = NULL
  
  return(output)
  
}
# wrapper for soma --------------------------------------------------------

.soma = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  control = NULL # check!
  
  con = soma::all2one()
  control = check_control(control=control, default=con)
  
  # we pass 'par' as initial population (gaussian with uniform-like variance).
  init = .createRandomPopulation(n=control$NP, mean=par, lower=lower, upper=upper) 
  
  xoutput = suppressWarnings(soma::soma(costFunction=fn, 
                                        bounds=list(min=lower, max=upper), 
                                        options=control, init=init))
  
  output = list()
  output$par  = xoutput$population[, xoutput$leader]
  output$value = xoutput$cost[xoutput$leader]
  output$counts = c('function'=xoutput$migrations*length(xoutput$cost),
                    gradient=NA)
  output$convergence = NA
  output$message = NULL
  
  return(output)
  
}


# wrapper for genoud ------------------------------------------------------

.genoud = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  print.level = if(control$verbose) control$trace else 0
  
  output = suppressWarnings(rgenoud::genoud(fn=fn, nvars=length(par), starting.values=par, 
                                            Domains = cbind(lower, upper),
                                            print.level=print.level))
  
  output$counts = c('function'=output$popsize*output$generations, gradient=NA)
  return(output)
  
}

# wrapper for pso ---------------------------------------------------------

.pso = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  hybrid = if(method=="hybridPSO") TRUE else FALSE
  method = if(any(method=="PSO", method=="PSO2007", method=="hybridPSO")) 
    "SPSO2007" else "SPSO2011"
  
  control = c(control, type=method, hybrid=hybrid)
  
  output = suppressWarnings(pso::psoptim(par=par, fn=fn, gr=gr, lower=lower, 
                                         upper=upper, control=control))
  
  return(output)
  
}


