.calibrar = function (par, fn, gr = NULL, ..., lower = -Inf, upper = Inf, active=NULL, 
                      control = list(), hessian = FALSE, method = NULL, skeleton=NULL,
                      replicates=1) {
  
  skeleton = skeleton
  if(is.null(skeleton)) skeleton = as.relistable(par)
  
  npar = length(par)
  
  parallel = control$parallel
  
  # check active parameters
  active = .checkActive(active=active, npar=npar)
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  # par can be re-scaled
  parscale = .checkParscale(control=control, npar=npar)
  control$parscale = NULL # reset parscale
  
  par   = par/parscale
  lower = lower/parscale
  upper = upper/parscale
  
  # function can be re-scaled
  fnscale = .checkFnscale(control=control)
  control$fnscale = NULL # reset fnscale
  
  # update to active parameters only
  guess  = par # the full original set of parameters
  par    = guess[isActive]
  lower  = lower[isActive]
  upper  = upper[isActive]
  
  npar = length(par)
  
  force(replicates)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  if(!is.null(gr)) gr = match.fun(gr)
  
  run = control$run
  
  # here we modify 'fn' so:
  # 0. it is evaluated in the '$run/..i' folder if necessary
  # 1. use 'isActive' to mask some parameters ('guess' is reference)
  # 2. is re-listed according to skeleton
  # 3. is re-scaled according to control$fnscale and control$parscale
  # 4. is evaluated a 'replicates' number of times
  fn1  = function(par, ..i=0) {
    
    if(!is.null(run)) {
      cwd = getwd()
      on.exit(setwd(cwd))
      .setWorkDir(run, i=..i) 
    }
    
    parx = guess
    parx[isActive] = par
    parx = relist(flesh = parx*parscale, skeleton = skeleton)
    output = NULL
    for(i in seq_len(replicates)) {
      out = fn(parx, ...)/fnscale
      output = rbind(output, t(as.matrix(out)))
    }
    if(is.null(output)) return(NULL)
    return(as.numeric(colMeans(output)))
    
  }
  
  attr(fn1, "..i") = TRUE
  
  # GRADIENT
  gr.control = control$gradient
  gr.method  = control$gr.method # if NULL, uses default.
  
  if(!is.null(gr)) {
    # here we modify 'gr' so:
    # 0. No modification, we assume gradient is computed in memory.
    # 1. use 'isActive' to mask some parameters ('guess' is reference)
    # 2. is re-listed according to skeleton
    # 3. is re-scaled according to control$fnscale and control$parscale
    # 4. No replication: if provided, we assume is deterministic.
    gr1  = function(par) {
      parx = guess
      parx[isActive] = par
      parx = relist(flesh = parx*parscale, skeleton = skeleton)
      grad = gr(parx, ...)*parscale/fnscale
      grad = grad[isActive]
      return(grad)
    }
  } else {
    gr1  = function(par) {
      gradient(fn=fn1, x=par, method=gr.method, control=gr.control, parallel=parallel)
    }    
  }
  
  # this need to be added so functions can use disk (specially in parallel)
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  n = .batchsize(par=par, method=method, control=control, parallel=parallel)
  copy_master_folder(control, n=n) # set a copy of master for all individuals
  
  if(all(method %in% multiMethods)) {
    # only run it if not provided in fn
    if(is.null(control$nvar)) control$nvar = length(fn1(par))
  }
  
  # here, make methods explicit (one by one)
  output = .all_methods(method=method, par=par, fn=fn1, gr=gr1, lower=lower, upper=upper, control=control, hessian=hessian)
  # reshaping full parameters, un-scaling par
  paropt = guess
  paropt[isActive] = output$par
  paropt = paropt*parscale
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")

  output$par = paropt[isActive]
  output$value = output$value*fnscale
  # here to rename 'par' to '.par'
  names(output)[names(output)=="par"] = ".par"
  # final outputs
  output = c(list(par=paropt), output, list(active=list(par=isActive, flag=activeFlag)))
  
  return(output)
  
}
