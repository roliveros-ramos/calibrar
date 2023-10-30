
.ahres = function(par, fn, gr, lower, upper, control, method, hessian) {
  
  # get restart for the current phase
  restart = .restartCalibration(control) # flag: TRUE or FALSE
  
  # default control options
  npar = length(par)
  con = list(trace = 0, fnscale = 1, parscale = rep.int(1L, npar), maxit=10000+200*npar,
             abstol = -Inf, reltol = sqrt(.Machine$double.eps), REPORT = 10L, 
             alpha=0.05, age.max=1, selection=0.5, step=0.5, weights=1,
             aggFn=.weighted.sum, nvar=1, useCV=TRUE, maxgen=NULL,
             convergence=1e-6, ncores=1, parallel=FALSE, verbose=FALSE,
             termination = 2, max_no_improvement=10, fn_smoothing=5)
  popOpt = .optPopSize(n=length(par), selection=con$selection)
  con$popsize = popOpt
  # end of default control options
  
  if(isTRUE(restart)) {
    
    res = .getRestart(control=control)
    opt   = res$opt
    trace = res$trace
    .messageByGen(opt, trace, restart=TRUE)
    
  } else {
    
    control = .check_control_ahres(control=control, default=con)
    
    opt = .newOpt(par=par, lower=lower, upper=upper, control=control)
    
    trace = NULL
      
    if(control$REPORT>0 & control$trace>0) {
      
      trace = list()
      trace$control = control
      trace$par   = matrix(NA, nrow=control$maxgen, ncol=length(par))
      trace$value = opt$the_values
      trace$best  = rep(NA, control$maxgen)
      
      trace$timing = list()
      trace$timing$fit   = rep(NA, control$maxgen)
      trace$timing$trace = rep(NA, control$maxgen)
      trace$timing$cont  = rep(NA, control$maxgen)
      
      if(control$termination %in% c(2,3)) {
        trace$stop  = rep(NA, control$maxgen)
        trace$timing$stop = rep(NA, control$maxgen)
      }
      
      if(control$trace>1) {
        trace$sd   = matrix(NA, nrow=control$maxgen, ncol=length(par))   
        trace$step = rep(NA, control$maxgen)     
      }
      
      if(control$trace>2) trace$fitness = NULL
      
      if(control$trace>3) trace$opt = vector("list", control$maxgen)
      
    }
  }
  
  # copy master folder after optimizing popsize
  copy_master_folder(control, n=control$popsize) 
  # start new optimization
  while(isTRUE(cc <- .continueEvolution(opt, control))) {
    
    tm1 = Sys.time()
    
    opt$gen  = opt$gen + 1
    opt$ages = opt$ages + 1
    
    # 1. create a new population
    if(all(opt$SIGMA==0)) break
    opt$pop = .createPopulation(opt)
    
    # 2. evaluate the function in the population: evaluate fn, aggregate fitness
    opt$fitness = .calculateFitness(opt, fn=fn)
    
    # 3. select best 'individuals'
    opt$selected = .selection(opt)
    
    # 4. create the new parents: MU and SD
    opt = .calculateOptimalWeights(opt)
    opt = .updatePopulation(opt)
    
    # save detailed outputs
    opt$the_values[opt$gen] = control$aggFn(opt$fitness[1, ], control$weights)
    
    xtm1 = Sys.time()
    opt$sstop[opt$gen] = smooth_stop(opt, control)
    xtm2 = Sys.time()
    tm_stop = format_difftime(xtm1, xtm2, value = TRUE)
    
    if(control$REPORT>0 & control$trace>0) {
      
      ttm1 = Sys.time()
      
      trace$par[opt$gen, ] = opt$MU
      trace$best[opt$gen]  = opt$selected$best$fit.global
      trace$value[opt$gen] = opt$the_values[opt$gen]
      trace$stop[opt$gen]  = opt$sstop[opt$gen]
      
      if(control$trace>1) {
        trace$sd[opt$gen, ]  = opt$SIGMA
        trace$step[opt$gen]  = opt$step       
      }
      
      if(control$trace>2) {
        
        if(is.null(trace$fitness)) 
          trace$fitness = matrix(NA, nrow=control$maxgen, ncol=ncol(opt$fitness))
        
        if(nrow(trace$fitness) < control$maxgen) {
          .nt = nrow(trace$fitness)
          trace$fitness = rbind(trace$fitness, 
                                matrix(NA, nrow=control$maxgen-.nt, ncol=ncol(opt$fitness))) 
        }
        
        trace$fitness[opt$gen, ] = opt$fitness[1, , drop=FALSE]
        
      }
      
      if(opt$gen %% control$REPORT==0) {
        if(control$trace>4) trace$opt[[opt$gen]] = opt
      }
      
      ttm2 = Sys.time()
      
      if(control$trace>3) {
        
      trace$timing$fit[opt$gen] = format_difftime(tm1, xtm1, value = TRUE)
      trace$timing$trace[opt$gen] = format_difftime(ttm1, ttm2, value = TRUE)
      trace$timing$cont[opt$gen] = max(attr(cc, "elapsed"),0) 
      trace$timing$stop[opt$gen] = tm_stop
      
      }
      
    }
    
    # save restart
    .createRestartFile(opt=opt, trace=trace, control=control)
    tm2 = Sys.time()
    
    if(control$verbose & opt$gen%%control$REPORT==0) {
      .messageByGen(opt, trace, level=control$trace, long=format_difftime(tm1, tm2))
    }
    
    
  } # end generations loop
  
  # trim trace
  last_gen = control$maxgen == opt$gen
  
  trace = .trim_trace_ahr(trace, n=opt$gen)
  
  use_disk = ifelse(!is.null(attr(fn, "..i")), TRUE, FALSE) 
  partial = if(use_disk) fn(opt$MU, ..i=0) else fn(opt$MU) 
  
  value = control$aggFn(x=partial, w=control$weights) # check if necessary
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, gradient=0)
  
  msg = sprintf("Stopping criteria reached in %d generations.", opt$gen)
  convergence = 0
  
  if(last_gen) {
    msg = "Maximum number of generations or function evaluations reached."
    convergence = 1
  }
  
  if(length(partial)==1) {

    output = list(par=opt$MU, value=value, counts=opt$counts, 
                  convergence=convergence, message=msg, generations=opt$gen)
    
  } else {

    output = list(par=opt$MU, value=value, counts=opt$counts, 
                  convergence=convergence, message=msg, generations=opt$gen, 
                  partial=partial)
    
  }
  
  if(!is.null(trace)) {
    output = c(output, trace=list(trace))
  }
    
  return(output)
  
}


# Auxiliar functions ------------------------------------------------------

# Initialize population ---------------------------------------------------

.calculateSigma = function(x, control=list()) {
  sigma = control$sigma
  if(is.null(control$sigma)) sigma = rep(NA, length(x))
  out = x^2/12
  out[!is.finite(out)] = 1
  sigma[is.na(sigma)] = out[is.na(sigma)]
  return(sigma)
}

.calculateRange = function(lower, upper) {
  out = upper - lower
  out[!is.finite(out)] = 1
  return(out)
}


.createRandomPopulation = function(n, mean, lower, upper) {
  mean = as.numeric(mean)
  lower = as.numeric(lower)
  upper = as.numeric(upper)
  SD	  = .calculateSigma(.calculateRange(lower, upper))
  out = rtnorm2(n=n, mean=mean, sd=SD, lower=lower, upper=upper)
  return(out)
}

.newOpt = function(par, lower, upper, control) {
  
  opt = list()
  opt$gen       = 0
  opt$npar      = length(par)
  opt$nvar      = control$nvar
  opt$seed      = control$popsize
  opt$selection = control$selection 
  opt$range     = .calculateRange(lower, upper)
  opt$MU        = as.numeric(par)
  opt$SIGMA	    = .calculateSigma(opt$range, control)
  opt$lower     = lower
  opt$upper     = upper
  opt$ps		    = numeric(opt$npar)
  opt$pc		    = numeric(opt$npar)
  opt$ages	    = numeric(control$popsize)
  opt$step	    = control$step   
  opt$mu		    = matrix(opt$MU, nrow=opt$npar, ncol=opt$nvar)
  opt$sigma	    = matrix(0, nrow=opt$npar, ncol=opt$nvar)
  
  opt$SD        = opt$step*sqrt(opt$SIGMA)
  
  opt$parents   = ceiling(opt$selection*opt$seed)
  opt$w.rec     = .getRecombinationWeights(parents=opt$parents, method=control$method)
  opt$mu.eff    = 1/sum(opt$w.rec*opt$w.rec)
  opt$cc        = 4/(opt$npar+4)
  opt$chiN      = sqrt(opt$npar)*(1-1/(4*opt$npar)+1/(21*opt$npar^2))
  
  opt$cs        = (opt$mu.eff+2)/(opt$npar+opt$mu.eff+3)
  opt$D         = 1 + 2*max(0, sqrt((opt$mu.eff-1)/(opt$npar+1))-1) + opt$cs
  opt$mu.cov    = opt$mu.eff
  opt$c.cov     = (1/opt$mu.cov)*(2/((opt$npar + sqrt(2))^2)) + 
    (1-1/opt$mu.cov)*min(1,(2*opt$mu.eff-1)/((opt$npar+2)^2+opt$mu.eff))
  
  opt$aggFn     = match.fun(control$aggFn)
  opt$weights   = control$weights
  opt$useCV     = control$useCV
  opt$alpha     = control$alpha
  
  opt$the_values = rep(NA, control$maxgen)
  opt$sstop      = rep(NA, control$maxgen)
  
  opt$control   = control
  
  class(opt) = c("restart", class(opt))
  return(opt)
  
}

.createPopulation = function(opt) {
  
  out = rtnorm2(n=opt$seed, mean=opt$MU, sd=opt$SD, lower=opt$lower, upper=opt$upper)
  out[, 1] = opt$MU # add the mean, it doesn't change the expected average of the population
  
  return(out)

}

.getRecombinationWeights = function(parents, method) {
  out = log(parents+1)-log(1:parents)
  out = out/sum(out)
  return(out)
}


# Calculate Fitness -------------------------------------------------------

  
.calculateFitness = function(opt, fn) {
  
  fn = match.fun(fn)
  pop = opt$pop
  gen = opt$gen
  parallel = opt$control$parallel
  run      = opt$control$run
  nvar     = opt$control$nvar
  
  use_disk = ifelse(!is.null(attr(fn, "..i")), TRUE, FALSE) 
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  if(isTRUE(parallel)) {
    # optimize parallel execution, reduce data transfer
    
    FITNESS  =  foreach(i=1:opt$seed, .verbose=FALSE, .inorder=FALSE) %dopar% {

      Fitness = if(use_disk) fn(pop[, i], ..i=i) else fn(pop[, i])  # internally set to ith wd
      if(is.null(Fitness)) {
        .write_calibrar_dump(run=run, gen=gen, i=i)
      }
        
      Fitness = c(i, Fitness)
      Fitness
      
    }
    
    FITNESS = .rbind_fitness(FITNESS)
    FITNESS = FITNESS[order(FITNESS[,1]), ][,-1, drop=FALSE]
    
  } else {
    
    FITNESS	=	NULL
    
    for(i in 1:opt$seed) {
      
      Fitness = if(use_disk) fn(pop[, i], ..i=i) else fn(pop[, i])
      if(is.null(Fitness)) {
        .write_calibrar_dump(run=run, gen=gen, i=i)
        Fitness = Inf
      }
      FITNESS	=	rbind(FITNESS, Fitness)
    }
    
  }
  
  return(FITNESS)
}


# Selection ---------------------------------------------------------------


.selection = function(opt) {
  
  pop     = opt$pop
  fitness = opt$fitness
  
  fitness.global = .globalFitness(fitness=fitness, opt=opt)
  
  p = sort(fitness.global,index.return=TRUE)
  
  supsG			= p$ix[seq_len(opt$parents)]
  .best     = p$ix[1]
  
  best = list(BEST       = pop[, .best],
              fit	       = fitness[.best, ],
              fit.global = fitness.global[.best])
  
  supsL = apply(fitness[supsG, ,drop=FALSE], 2, FUN = order)
  
  return(list(supsG=supsG, supsL=supsL, best=best))
  
}


# Recombination -----------------------------------------------------------

.norma = function(x) {
  # x[x==0] = 1E-20
  return(x/sum(x, na.rm=TRUE))
  
}

.w.oi=function(x, b=4) {
  
  n       =	ncol(x)
  if(n==1) return(matrix(1, nrow=nrow(x), ncol=1))
  
  cv.min  = 0.9*(apply(x, 2, min, na.rm=TRUE) + 1e-20)
  cv.max  = 1.1*(apply(x, 2, max, na.rm=TRUE) + 1e-20)
  out     = t(x)
  out     = ((cv.max - out)/(cv.max-cv.min))^b
  out     = t(out/rowSums(out))
  out[out==0] = 1e-20
  out	    = t(apply(out, 1, .norma))
  return(out)
  
}


.calculateOptimalWeights = function(opt) {
  
  pop   = opt$pop[, opt$selected$supsG]
  supsL = opt$selected$supsL
  
  w.rec = opt$w.rec
  alpha = opt$alpha
  
  opt.ind	= array(NA, dim=c(nrow(pop), ncol(supsL)))
  opt.sd	= array(NA, dim=c(nrow(pop), ncol(supsL)))
  
  for(i in seq_len(ncol(supsL))) {
    pop_i          = t(pop[, supsL[, i]])*w.rec
    opt.ind[, i]   = colSums(pop_i)
    opt.var	       = colSums(pop_i^2) - opt.ind[, i]^2
    # pop_i          = pop[, supsL[, i]])
    # opt.ind[, i]   = apply(pop_i, 1, weighted.mean, w=w.rec)
    # opt.var	       = apply(pop_i^2, 1, weighted.mean, w=w.rec) - opt.ind[, i]^2
    opt.var        = pmax(opt.var, 0)
    opt.sd[, i]    = sqrt(opt.var)
  }
  
  mu.new		= (1-alpha)*opt$mu + alpha*opt.ind
  s.new			= (1-alpha)*(opt$mu^2+opt$sigma^2) + alpha*(opt.ind^2+ opt.sd^2) - mu.new^2
  s.new     = pmax(s.new, 0) # s.new is always positive, this is a correction for rounding errors 
  sigma.new = sqrt(s.new)
  ww.new		= if(isTRUE(opt$useCV)) .w.oi(sigma.new/opt$range) else .w.oi(sigma.new)
  
  ww.rec = array(w.rec[supsL], dim=dim(supsL))
  W      = ww.new %*% t(ww.rec)

  opt$W     = W
  opt$w.new = ww.new
  opt$mu    = mu.new
  opt$sigma = sigma.new
  opt$best  = opt$selected$best # updateBestIndividual()
  
  opt$MU.eff  = 1/rowSums(W*W)
  opt$mu.eff  = max(mean(opt$MU.eff), 1) # mu.eff is always greater or equal to 1, avoiding rounding errors
    
  return(opt)
  
}


.updatePopulation = function(opt) {
  
    supsG       = opt$selected$supsG
    npar        = opt$npar
    MU          = opt$MU
    SIGMA       = opt$SIGMA
    step        = opt$step
    
    opt$cs      = (opt$mu.eff+2)/(npar+opt$mu.eff+3)
    opt$D       = 1 + 2*max(0, sqrt((opt$mu.eff-1)/(npar+1))-1) + opt$cs
    opt$mu.cov  = opt$mu.eff
    opt$c.cov   = (1/opt$mu.cov)*(2/((npar+sqrt(2))^2))+(1-1/opt$mu.cov)*
      min(1,(2*opt$mu.eff-1)/((npar+2)^2+opt$mu.eff))
    
    MU.new      = rowSums(opt$W * opt$pop[, supsG, drop=FALSE])
    
    opt$pc      = (1-opt$cc)*opt$pc + sqrt(opt$cc*(2-opt$cc))*sqrt(opt$MU.eff)*(MU.new-MU)/step
    opt$ps      = (1-opt$cs)*opt$ps + sqrt(opt$cs*(2-opt$cs))*sqrt(opt$MU.eff)*((MU.new-MU)/sqrt(SIGMA))/step
    
    SIGMA.sel   = rowSums(opt$W * ((opt$pop[, supsG, drop=FALSE] - MU)/step)^2)
    SIGMA.new   = (1-opt$c.cov)*SIGMA + (opt$c.cov/opt$mu.cov)*opt$pc*opt$pc + 
      opt$c.cov*(1-1/opt$mu.cov)*SIGMA.sel
    
    opt$step    = step*exp(opt$cs*(sqrt(sum(opt$ps*opt$ps, na.rm=TRUE))/opt$chiN-1)/opt$D)
    
    opt$MU      = MU.new
    opt$SIGMA   = SIGMA.new
    
    opt$SD      = opt$step*sqrt(opt$SIGMA)
    
  return(opt)
  
}


# Auxiliar functions ------------------------------------------------------

.check_control_ahres = function(control, default) {
  
  popOpt = default$popsize
  control = check_control(control=control, default=default, minimal=FALSE, verbose=FALSE)
  # update population size and selection rate
  if(control$popsize < popOpt) {
    warning("'popsize' is too small, using default value.")
    message(sprintf("Population size set to %d (optimal value).", popOpt))
  }
  control$popsize = max(control$popsize, popOpt)
  
  if(isTRUE(control$parallel)) {
    popOptP = ceiling(control$popsize/control$ncores)*control$ncores
    if(control$popsize != popOptP) {
      message(sprintf("Optimizing 'popsize' to work with %d cores...", control$ncores))
      message(sprintf("\tOptimal population size (%d) has been increased to %d.\n", popOpt, popOptP))
    }
    
    control$popsize = popOptP
    control$selection = control$selection*popOpt/popOptP
  }
  # update maximum number of function evaluations and generations
  if(!is.null(control$maxit) & !is.null(control$maxgen))
    message("'maxit' and 'maxgen' provided, ignoring 'maxgen'.")
  
  if(is.null(control$maxit) & is.null(control$maxgen))
    stop("You must provide 'maxit' or 'maxgen' for 'AHR-ES' method.")
  
  if(is.null(control$maxit) & !is.null(control$maxgen))
    control$maxit = control$maxgen
  
  control$maxgen = control$maxit
  
  return(control)
  
}


.trim_trace_ahr = function(trace, n) {
  if(is.null(trace)) return(NULL)
  ind = seq_len(n)
  trace$par     = trace$par[ind, ]
  trace$value   = trace$value[ind]
  trace$best    = trace$best[ind]
  trace$stop    = trace$stop[ind]
  trace$sd      = trace$sd[ind, ]
  trace$step    = trace$step[ind]
  trace$fitness = trace$fitness[ind, ]
  trace$opt     = trace$opt[ind]
  if(!is.null(trace$timing)) {
    trace$timing$fit   = trace$timing$fit[ind]
    trace$timing$trace = trace$timing$trace[ind]
    trace$timing$cont  = trace$timing$cont[ind]
    trace$timing$stop  = trace$timing$stop[ind]
  }
  return(trace)
}
