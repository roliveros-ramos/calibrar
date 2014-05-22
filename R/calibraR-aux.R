
.checkActive = function(active, npar) {

  if(is.null(active)) return(rep(TRUE, npar))
  if(all(!active)) stop("No parameter is active, at least one parameter must be optimized.")
  if(!is.logical(active)) stop("'active' must be a boolean vector.")
  if(length(active)!=npar) stop("Length of 'active' parameters must match parameter number.")
  active[is.na(active)] = FALSE
  return(active)
  
}

.checkBounds =  function(lower, upper, npar) {
  
  nl = length(lower)
  nu = length(upper)

  if(nl==1) {
    if(npar!=1 & lower!=-Inf) 
      warning("Only one lower bound has been provided, used for all parameters.")
    lower = rep(lower, npar)
    nl = npar
  }
  
  if(nu==1) {
    if(npar!=1 & upper!=Inf) 
      warning("Only one upper bound has been provided, used for all parameters.")    
    upper = rep(upper, npar)
    nu = npar
  }
  
  if(nl!=npar) stop("Lower bounds must match parameter number.")
  if(nu!=npar) stop("Upper bounds must match parameter number.")
  
  if(any(is.na(lower))) {
    lower[is.na(lower)] = -Inf
    warning("NAs supplied in lower thresholds, replaced by -Inf")
  }
  
  if(any(is.na(upper))) {
    upper[is.na(upper)] = Inf
    warning("NAs supplied in upper thresholds, replaced by Inf")
  }
  
  if(any(lower >= upper)) stop("Lower bounds must be lower than upper bounds.")
  
  output = list(lower=lower, upper=upper)
  
  return(output)
  
}

.checkOpt = function(par, lower, upper) {

  if(any(par < lower, na.rm=TRUE)) 
    stop("Initial guess for parameters must be greater than lower threshold.")
  
  if(any(par > upper, na.rm=TRUE)) 
    stop("Initial guess for parameters must be lower than upper threshold.")
  
  ind = is.na(par) & !is.finite(lower) & !is.finite(upper)
  par[ind] = 0
  
  ind = is.na(par) & is.finite(lower) & is.finite(upper)
  par[ind] = 0.5*(lower + upper)[ind]
  
  ind = is.na(par) & is.finite(lower) & !is.finite(upper)
  par[ind] = lower[ind] + 0.1*abs(lower[ind])
  
  ind = is.na(par) & !is.finite(lower) & is.finite(upper)
  par[ind] = upper[ind] - 0.1*abs(upper[ind]) 
  
  return(par)
  
}

.checkControl = function(control, method, par, fn, active) {
  
  fn = match.fun(fn)
  
  con = list(trace = 0, fnscale = 1, parscale = rep.int(1, length(par)), maxit = NULL, maxgen=NULL,
             abstol = -Inf, reltol = sqrt(.Machine$double.eps), REPORT = 5, ncores=parallel::detectCores(), 
             alpha=0.05, age.max=1, selection=0.5, step=0.5, nvar=NULL, weights=1, sigma=NULL,
             method=method, aggFn=.weighted.sum, parallel=FALSE, run=NULL, useCV=TRUE,
             convergence=1e-6)
  
  popOpt = floor(0.5*(4 + floor(3*log(length(par))))/con$selection)
  con$popsize = popOpt
  
  controlDef  = names(con)       # default options
  controlUser = names(control)   # user provided options
  
  con[controlUser] = control
  
  # get unknown options
  unknown = controlUser[!(controlUser %in% controlDef)]
  
  # update population size and selection rate
  if(con$popsize < popOpt) warning("'popsize' is too small, using default value.")
  con$popsize = max(con$popsize, popOpt)
  if(isTRUE(con$parallel)) {
    popOptP = ceiling(con$popsize/con$ncores)*con$ncores
    if(con$popsize != popOptP) message(paste("Optimizing 'popsize' to work with", con$ncore, "cores."))
    con$popsize = popOptP
    con$selection = round(con$selection*popOpt/popOptP, 1)
  }
  
  # check for user-provided variances (sigma)
  if(!is.null(con$sigma) & length(con$sigma)!=length(active)) 
    stop("Vector of variances (sigma) must match parameter length.")
  if(!is.null(con$sigma)) con$sigma = con$sigma[which(active)]
  
  # check number of variables
  if(is.null(con$nvar)) con$nvar = length(fn(par))
  
  # update maximum number of function evaluations and generations
  if(!is.null(con$maxit) & !is.null(con$maxgen)) 
    warning("'maxit' and 'maxgen' provided, ignoring 'maxit'.")
  if(!is.null(con$maxit) & is.null(con$maxgen)) con$maxgen = floor(con$maxit/con$popsize)
  if(is.null(con$maxit) & is.null(con$maxgen)) con$maxgen = 1000L
  
  con$maxit = con$popsize*con$maxgen
  
  # update and check weights
  if(length(con$weights)==1) con$weights = rep(con$weights, con$nvar)
  if(length(con$weights)!=con$nvar) stop("Vector of weights should match the length of the output of fn.")
  if(any(con$weights<0)) stop("Weights should be positive numbers.")
  if(any(is.na(con$weights))) stop("Weights cannot be NA.")
  
  # aggregation function for global fitness
  con$aggFn = match.fun(con$aggFn)
  
  if(length(unknown)!=0) 
    warning("Unknown control parameters: ", paste0(paste(unknown, collapse = ", ")),".")
  
  return(con)

}

.setWorkDir = function(run, i) {
  
  if(is.null(run)) return(invisible())
  work.dir = file.path(run, paste0("i", i))
  setwd(work.dir)
  return(work.dir)
  
}

.getRestart = function(restart) {
  
  return(restart)
}

.createRestartFile = function() {
  # MU, SIGMA, sigmasmall, gen, phase
}


