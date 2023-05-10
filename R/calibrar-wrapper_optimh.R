
# wrapper for cma_es ------------------------------------------------------

.cmaes = function(par, fn, lower, upper, control) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  npar = length(unlist(par))
  output = suppressWarnings(cmaes::cma_es(par=par, fn=fn, lower=lower, upper=upper, control=control))
  
  if(is.null(output$par)) {
    output$par = relist(rep(NA, npar), skeleton=par)
    if(!is.finite(output$value)) warning("Infinite value reached for fn.")
  }
  
  # names(output)[names(output)=="par"] = ".par"
  
  return(output)
  
}



# wrapper for genSA -------------------------------------------------------

.genSA = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  output = suppressWarnings(GenSA::GenSA(par=par, fn=fn, lower=lower, 
                                         upper=upper, control=control))
  
  # names(output)[names(output)=="par"] = ".par"
  output$counts = c('function'=output$counts,
                    gradient=NA)
  
  return(output)
  
}
# wrapper for DEoptim -----------------------------------------------------

.DE = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  if(isTRUE(control$parallel)) control$parallelType=2
  
  ctrl = DEoptim::DEoptim.control()
  control[!(names(control) %in% names(ctrl))] = NULL
  
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
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  control = NULL # check!
  
  xoutput = suppressWarnings(soma::soma(costFunction=fn, 
                                        bounds=list(min=lower, max=upper), 
                                        options=control))
  
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
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  print.level = if(control$verbose) control$trace else 0
  
  output = suppressWarnings(rgenoud::genoud(fn=fn, nvars=length(par), starting.values=par, 
                                            Domains = cbind(lower, upper),
                                            print.level=print.level))
  
  # names(output)[names(output)=="par"] = ".par"
  output$counts = c('function'=output$popsize*output$generations, gradient=NA)
  return(output)
  
}

# wrapper for pso ---------------------------------------------------------

.pso = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  hybrid = if(method=="hybridPSO") TRUE else FALSE
  method = if(any(method=="PSO", method=="PSO2007", method=="hybridPSO")) 
    "SPSO2007" else "SPSO2011"
  
  control = c(control, type=method, hybrid=hybrid)
  
  output = suppressWarnings(pso::psoptim(par=par, fn=fn, gr=gr, lower=lower, 
                                         upper=upper, control=control))
  
  # names(output)[names(output)=="par"] = ".par"
  
  return(output)
  
}


