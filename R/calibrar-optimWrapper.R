.calibrar = function (par, fn, gr = NULL, ..., lower = -Inf, upper = Inf, active=NULL, 
                      control = list(), hessian = FALSE, method = NULL, skeleton=NULL,
                      replicates=1) {
  
  skeleton = skeleton
  if(is.null(skeleton)) skeleton = as.relistable(par)
  
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
  
  force(replicates)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  if(!is.null(gr)) gr = match.fun(gr)
  
  control = .checkControl(control=control, method=method, par=guess, fn=fn, 
                          active=active, skeleton=skeleton, 
                          replicates=replicates, ...)
 
  if(is.null(control$master)) {
    
    fn1  = function(par, .i=0) {
      
      parx = guess
      parx[isActive] = par
      parx = relist(flesh = parx, skeleton = skeleton)
      output = NULL
      for(i in seq_len(replicates)) {
        out = fn(parx, ...)/control$fnscale
        output = rbind(output, t(as.matrix(out)))
      }
      return(as.numeric(colMeans(output)))
      
    }
    
  } else {

    fn1  = function(par, .i=0) {
      
      cwd = getwd()
      on.exit(setwd(cwd))
      .setWorkDir(control$run, i=.i) 
      
      parx = guess
      parx[isActive] = par
      parx = relist(flesh = parx, skeleton = skeleton)
      output = NULL
      for(i in seq_len(replicates)) {
        out = fn(parx, ...)/control$fnscale
        output = rbind(output, t(as.matrix(out)))
      }
      return(as.numeric(colMeans(output)))
      
    }
    
  } # end fn1
  
  
  imethod = "default"
  optimMethods  = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                    "Brent")
  optimxMethods = c("nlm", "nlminb", "spg", "ucminf", "newuoa", "bobyqa", 
                    "nmkb", "hjkb")
  psoMethods    = c("PSO", "PSO2007", "PSO2011", "hybridPSO")
  
  if(method %in% optimMethods) imethod = "optim"
  if(method %in% optimxMethods) imethod = "optimx"
  if(method %in% psoMethods) imethod = "pso"
  if(method == "Rcgmin") imethod = "Rcgmin"  
  if(method == "Rvmmin") imethod = "Rvmmin"  
  if(method == "cmaes") imethod = "cmaes"  
  if(method == "lbfgsb3") imethod = "lbfgsb3"  
  if(method == "genSA") imethod = "genSA"
  if(method == "soma") imethod = "soma"
  if(method == "genoud") imethod = "genoud"
  if(method == "DE") imethod = "DEoptim"
  
  output = 
    switch(imethod, 
           default = .optimES(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                              hessian=hessian, method=method, isActive=isActive),
           optim   = .optim(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                            hessian=hessian, method=method),
           optimx  = .optimx(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                             hessian=hessian, method=method),
           Rcgmin  = .Rcgmin(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                             hessian=hessian, method=method),
           Rvmmin  = .Rvmmin(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                             hessian=hessian, method=method),
           pso     = .pso(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                          hessian=hessian, method=method),
           cmaes   = .cmaes(par=par, fn=fn1, lower=lower, upper=upper, control=control),
           lbfgsb3 = .lbfgsb3(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                              hessian=hessian, method=method),
           genSA   = .genSA(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                            hessian=hessian, method=method),
           soma    = .soma(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                           hessian=hessian, method=method),
           genoud  = .genoud(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                             hessian=hessian, method=method),
           DEoptim = .DE(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                         hessian=hessian, method=method)
    )
  
  # reshaping full parameters
  paropt = guess
  paropt[isActive] = output$ppar 
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  # final outputs
  output = c(list(par=paropt), output, list(active=list(par=isActive, flag=activeFlag)))
  
  class(output) = c("optimES.result", class(output)) # change name
  
  return(output)
  
}


# wrapper for optim -------------------------------------------------------

.optim = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  # what to do with methods not accepting bound conditions?
  if(!(method %in% c("L-BFGS-B", "Brent"))) {
    lower = -Inf
    upper = Inf
  }
  
  output = suppressWarnings(stats::optim(par=par, fn=fn, gr=gr, method=method, lower=lower, 
                                         upper=upper, control=control, hessian=hessian))
  
  names(output)[names(output)=="par"] = "ppar"
  
  return(output)
  
}


# wrapper for optimx ------------------------------------------------------

.optimx = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  parNames = names(par)
  control = NULL # check!
  
  out = suppressWarnings(optimr::optimr(par=par, fn=fn, gr=gr, method=method, lower=lower, 
                                        upper=upper, control=control, hessian=hessian))
  
  par = as.numeric(out[, parNames])
  value = out$value
  counts = as.numeric(out[, c("fevals", "gevals")])
  output = list(ppar=par, value=value, counts=counts, convergence=out$convcode)
  
  return(output)
  
}


# wrapper for Rcgmin ------------------------------------------------------

.Rcgmin = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  parNames = names(par)
  control = NULL # check!
  
  output = suppressWarnings(Rcgmin::Rcgmin(par=par, fn=fn, gr=gr, lower=lower, 
                                           upper=upper, control=control))
  
  names(output)[names(output)=="par"] = "ppar"
  
  return(output)
  
}

# # wrapper for Rvmmin ----------------------------------------------------

.Rvmmin = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  parNames = names(par)
  control = NULL # check!
  
  output = suppressWarnings(Rvmmin::Rvmmin(par=par, fn=fn, gr=gr, lower=lower, 
                                           upper=upper, control=control))
  
  names(output)[names(output)=="par"] = "ppar"
  
  return(output)
  
}

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
  
  names(output)[names(output)=="par"] = "ppar"
  
  return(output)
  
}

# wrapper for lbfgsb3 -----------------------------------------------------

.lbfgsb3 = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  ctrl = list(maxit = 500, trace = 0, iprint = 0L)
  control[!(names(control) %in% names(ctrl))] = NULL
  
  xoutput = suppressWarnings(lbfgsb3::lbfgsb3(prm=par, fn=fn, gr=gr, lower=lower,
                                              upper=upper, control=control))
  
  output = list()
  output$ppar  = xoutput$prm
  output$value = xoutput$f
  output$counts = c('function'=xoutput$info$isave[34],
                    gradient=NA)# check on split!
  output$convergence = NA
  output$message = NULL
  
  return(output)
  
}

# wrapper for genSA -------------------------------------------------------

.genSA = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  output = suppressWarnings(GenSA::GenSA(par=par, fn=fn, lower=lower, 
                                         upper=upper, control=control))
  
  names(output)[names(output)=="par"] = "ppar"
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
  output$ppar  = xoutput$optim$bestmem
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
  output$ppar  = xoutput$population[, xoutput$leader]
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
  
  names(output)[names(output)=="par"] = "ppar"
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
  
  names(output)[names(output)=="par"] = "ppar"
  
  return(output)
  
}

# optimES internal --------------------------------------------------------

.optimES = function(par, fn, gr, lower, upper, control, method, hessian, isActive) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
 
  copyMaster(control) # set a copy of master for all individuals

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
      trace$par   = matrix(NA, nrow=control$maxgen, ncol=length(isActive))
      trace$value = rep(NA, control$maxgen)
      trace$best  = rep(NA, control$maxgen)
      
      if(control$trace>1) {
        trace$sd   = matrix(NA, nrow=control$maxgen, ncol=length(isActive))   
        trace$step = rep(NA, control$maxgen)     
      }
      
      if(control$trace>2) trace$fitness = NULL
      
      if(control$trace>3) trace$opt = vector("list", control$maxgen)
      
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
      
      if(control$trace>2) {
        
        if(is.null(trace$fitness)) 
          trace$fitness = matrix(NA, nrow=control$maxgen, ncol=ncol(opt$fitness)) 
        
        trace$fitness[opt$gen, ] = colMeans(opt$fitness[opt$selected$supsG, , drop=FALSE])
        
      }
      
      if(opt$gen%%control$REPORT==0) {
        if(control$trace>3) trace$value[opt$gen] = control$aggFn(fn(opt$MU), control$weights)
        if(control$trace>4) trace$opt[[opt$gen]] = opt
      }

    }
    
    # save restart
    .createRestartFile(opt=opt, trace=trace, control=control)
    
    if(control$verbose & opt$gen%%control$REPORT==0) 
      .messageByGen(opt, trace)
    
  } # end generations loop
  
  partial = fn(opt$MU, .i=0)
  value = control$aggFn(x=partial, w=control$weights) # check if necessary
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)
  
  output = list(ppar=opt$MU, value=value, counts=opt$counts, 
                trace=trace, partial=partial, convergence=1)
  
  
  return(output)
  
}


# Auxiliar functions to run fn on disk ------------------------------------

copyMaster = function(control, n=NULL) {
  if(is.null(n)) n = control$popsize
  for(i in (seq_len(n) - 1)) .copyMaster(control, i)
  return(invisible())
}

.copyMaster = function(control, i) {
  if(is.null(control$master)) return(invisible())
  if(is.null(control$run))
    stop("You must specify a 'run' directory to copy the contents of the 'master' folder.")
  newDir = .getWorkDir(control$run, i)
  if(!dir.exists(newDir)) dir.create(newDir, recursive=TRUE)
  if(isTRUE(control$master)) return(invisible(newDir))
  if(!dir.exists(control$master)) stop("The 'master' directory does not exist.") 
  xfiles = dir(path=control$master)
  from   = file.path(control$master, xfiles)
  file.copy(from=from, to=newDir, recursive=TRUE, overwrite = FALSE) # why?
  return(invisible(newDir))
}

.getWorkDir = function(run, i) {
  if(is.null(run)) return(invisible())
  work.dir = file.path(run, paste0("i", i))
  return(work.dir)
}
