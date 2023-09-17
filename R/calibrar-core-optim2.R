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
      gradient(fn=fn1, x=par, method=gr.method, control=gr.control, parallel=parallel)
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
