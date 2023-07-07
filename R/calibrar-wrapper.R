# all methods

.all_methods = function(method, par, fn, gr, lower, upper, control, hessian) {
  output = 
    switch(method, 
           "Nelder-Mead" = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="Nelder-Mead"), 
           "BFGS"        = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="BFGS"), 
           "CG"          = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="CG"), 
           "L-BFGS-B"    = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="L-BFGS-B"),
           "SANN"        = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="SANN"),
           "Brent"       = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="Brent"), 
           "nlm"         = .nlm(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "nlminb"      = .nlminb(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "Rcgmin"      = .Rcgmin(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "Rvmmin"      = .Rvmmin(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "hjn"         = .hjn(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "spg"         = .spg(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "LBFGSB3"     = .lbfgsb3(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "AHR-ES"      = .ahres(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           stop(printf("UNSUPPORTED METHOD: %s.", sQuote(method)), call. = FALSE)
    )
  return(output)
}

# calibrate ---------------------------------------------------------------

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
      gradient(fn=fn1, x=par, method=gr.method, control=gr.control, parallel=parallel, ...)
    }    
  }
  
  # this need to be added so functions can use disk (specially in parallel)
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  n = .batchsize(par=par, method=method, control=control)
  copy_master_folder(control, n=n) # set a copy of master for all individuals
  
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


# .optim2 wrapper ---------------------------------------------------------

.optim2 = function(par, fn, gr = NULL, ..., lower = -Inf, upper = +Inf, active = NULL, 
                   method = NULL, control = list(), hessian = FALSE, parallel=FALSE) {
  
  skeleton = as.relistable(par)
  
  par    = unlist(par)
  lower  = unlist(lower)
  upper  = unlist(upper)
  active = unlist(active)
  
  npar = length(par)
  
  # par can be re-scaled
  parscale = .checkParscale(control=control, npar=npar)
  control$parscale = NULL # reset parscale
  
  par   = par/parscale
  lower = lower/parscale
  upper = upper/parscale
  
  # function can be re-scaled
  fnscale = .checkFnscale(control=control)
  control$fnscale = NULL # reset fnscale
  
  # Checking active par and bounds
  active     = .checkActive(active=active, npar=npar)
  bounds     = .checkBounds(lower=lower, upper=upper, npar=npar)
  # update to active parameters only
  guess   = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  par     = guess
  lower   = bounds$lower
  upper   = bounds$upper
  
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  if(is.null(names(par))) names(par) = .printSeq(npar, preffix="par")
  
  # closure for function evaluation
  fn = match.fun(fn)
  if(!is.null(gr)) gr = match.fun(gr)
  
  # here we modify 'fn' so:
  # 1. use 'isActive' to mask some parameters ('guess' is reference)
  # 2. is re-listed according to skeleton
  # 3. is re-scaled according to control$fnscale and control$parscale
  
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    parx = relist(flesh = parx*parscale, skeleton = skeleton)
    output = fn(parx, ...)/fnscale
    return(output)
  }
  
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
      gradient(fn=fn1, x=par, method=gr.method, control=gr.control, parallel=parallel, ...)
    }    
  }
  
  # here, make methods explicit (one by one)
  output = .all_methods(method=method, par=par, fn=fn1, gr=gr1, lower=lower, upper=upper, control=control, hessian=hessian)
  
  # reshaping full parameters, un-scaling par
  paropt = guess
  paropt[isActive] = output$par
  paropt = paropt*parscale
  
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  paropt = relist(paropt, skeleton)
  class(paropt) = setdiff(class(paropt), "relistable")
  
  output$par = paropt
  output$value = output$value*fnscale
  
  return(output) 
  
}


# Manage control options --------------------------------------------------

check_control = function(control, default, minimal=TRUE, verbose=TRUE) {
  
  control = control[!sapply(control, is.null)] # remove NULL values
  nm_full = names(default)
  ignored = setdiff(names(control), nm_full)
  keep = names(default)[nm_full %in% names(control)]
  if(isTRUE(minimal)) return(control[keep])
  msg = sprintf("Ignoring control arguments: %s.", paste(sQuote(ignored), collapse=", "))
  if(length(ignored) & verbose) message(msg)
  default[keep] = control[keep]
  return(default)
}


# Auxiliar functions to run fn on disk ------------------------------------

copy_master_folder = function(control, n=NULL) {
  if(is.null(control$master)) return(invisible())
  if(is.null(control$run))
    stop("You must specify a 'run' directory to copy the contents of the 'master' folder.")
  if(is.null(n)) return(invisible())
  for(i in (seq_len(n+1) - 1)) .copyMaster(control, i)
  return(invisible())
}

.copyMaster = function(control, i) {
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
