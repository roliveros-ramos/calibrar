.setNewCalibration = function() {
  return(invisible())
}
.calculateObjetiveValue = function(obs, sim, info) {
  fit = NULL
  for(j in seq_len(nrow(info))) {
    if(!info$calibrate[j]) next
    var = info$variable[j]
    msg = "Variable '%s' not found in %s data (existing variables: %s)."
    if(!(var %in% names(obs))) 
      stop(sprintf(msg, var, "observed", paste(sQuote(names(sim)), collapse=", ")))
    if(!(var %in% names(sim)))
      stop(sprintf(msg, var, "simulated", paste(sQuote(names(sim)), collapse=", ")))

    fit = c(fit, .fitness(obs=obs[[var]], sim=sim[[var]], FUN=info$type[j]))
  }
  names(fit) = info$variable[which(info$calibrate)]
  return(fit)
}

.fitness = function(obs, sim, FUN, ...) {
  FUN = match.fun(FUN)
  output = FUN(obs=obs, sim=sim, ...)
  return(output)
}

.messageByGen = function(opt, trace, level=0, restart=FALSE, long=NULL) {
  
  nx  = ceiling(log10(opt$control$maxgen + 0.5))

  msg0 = paste0("Generation %", nx, "d finished (%s).\n")
  msg1 = paste0("Generation %", nx, "d finished (%s)\n   Function value = %.10g\n")
  msg2 = paste0("Generation %", nx, "d finished (from restart).\n")

  if(is.null(long)) long = date()
  
  if(isTRUE(restart)) {
    out = sprintf(msg2, opt$gen)
    message(out)
    return(invisible(NULL))
  }
  
  if(level>3) {
    out = sprintf(msg1, opt$gen, long, trace$value[opt$gen])
    message(out)
    return(invisible(NULL))
  }
  
  if(level>0) {
    out = sprintf(msg1, opt$gen, long, trace$best[opt$gen])
    message(out)
    return(invisible(NULL))
  }
  
  out = sprintf(msg0, opt$gen, long)
  message(out)
  
  return(invisible(NULL))
}

# .messageByGen = function(opt, trace) {
#   cat("--- Generation", opt$gen, "-----")
# }


.printSeq = function(n, preffix=NULL, suffix=NULL, sep="") {
  nx  = ceiling(log10(n))
  fmt = paste0(sprintf(paste0("%0", nx, "d"), seq_len(n)))
  out = paste(preffix, fmt, suffix, sep=sep)
  return(out)
}

.read.csv3 = function(path,  ...) {
  
  x = tryCatch(read.csv(file=path, row.names=1, ...),
               error = function(e) .errorReadingCsv(e, path))
  if(is.null(x)) stop()
  x = as.matrix(x)
  mode(x)="numeric"
  
  return(x)
}

.read_data  = function(file, col_skip=NA, varid=NA, nrows=-1, ...) {

  type = .guess_type(file)
  
  if(type=="text") {
    out = .read.csv4(file=file, col_skip=col_skip, varid=varid, nrows=nrows, ...)
    return(out)
  }
  
  if(type=="rds") {
    out = readRDS(file=file)
    return(out)
  }
  
  if(type=="nc") {
    message("netCDF files are not yet supported.")
  }
  
  stop(sprintf("File type not supported (%s)", file))  
  
}

.guess_type = function(file) {
  ext = tail((unlist(strsplit(file, split="\\."))), 1)
  is_txt = grepl(ext, pattern="csv|txt")
  is_nc = grepl(ext, pattern="nc")
  is_rds = grepl(ext, pattern="rds")
  if(is_txt) return("text")
  if(is_nc) return("nc")
  if(is_rds) return("rds")
}

.read.csv4 = function(file, col_skip=NA, varid=NA, nrows=-1, ...) {
  if(is.na(nrows)) nrows=-1
  x = tryCatch(read.csv(file=file, row.names=NULL, check.names = FALSE, nrows=nrows, ...),
               error = function(e) .errorReadingCsv(e, file))
  if(is.null(x)) stop(sprintf("File %s not found.", file))
  if(!is.na(varid)) {
    if(length(varid)!=1) stop("Only one varid must be provided.")
    if(!(varid %in% colnames(x))) 
      stop(sprintf("varid '%s' not found in %s", varid, file))
    return(as.numeric(x[, varid]))
  }
  x = as.matrix(x)
  if(!is.na(col_skip)) {
    x = x[, -seq_len(col_skip)]
  }
  mode(x)="numeric"
  return(x)
}


.errorReadingCsv = function(e, path) {
  on.exit(stop())
  message(e, ".")
  message("Error reading file ", sQuote(path), ".", sep="")
  return(invisible())
}


.pasteList = function(x) {
  if(length(x)<=1) return(x)
  out = paste(paste(x[-length(x)], sep="", collapse=", "), "and", x[length(x)])
  return(out)
  
}

format_difftime = function(x, y, ...) {
  dt = as.numeric(difftime(y, x, units="secs"))
  h = trunc(dt/3600)
  m = trunc((dt - h*3600)/60)
  s = dt - h*3600 - m*60
  
  if(h==0) {
    out = sprintf("%dm %0.1fs", m, round(s, 1))
    if(m==0) {
      out = sprintf("%0.2fs", round(s, 2))
    }
    return(out)
  }
  out = sprintf("%dh %dm %0.1fs", h, m, round(s, 1))
  return(out)
}

.checkActive = function(active, npar) {
  
  if(is.null(active)) return(rep(TRUE, npar))
  active = as.logical(active)
  if(all(!active)) stop("No parameter is active, at least one parameter must be optimized.")
  if(length(active)!=npar) stop("Length of 'phases' argument must match number of parameters.")
  return(active)
  
}

.checkPhases = function(phases, npar) {
  
  if(is.null(phases)) return(rep(1L, npar))

  phases = as.integer(phases)
  active = phases>0 & !is.na(phases)
  if(all(!active)) stop("No parameter is active, at least one parameter must be optimized.")
  if(length(active)!=npar) stop("Length of 'phases' argument must match number of parameters.")
  phases[phases<=0] = NA
  if(min(phases, na.rm=TRUE)!=1L) stop("No parameters has been specified for phase 1.")
  allPhases = seq(from=1, to=max(phases, na.rm=TRUE))
  ind = allPhases %in% phases
  missingPhases = .pasteList(allPhases[!ind])
  if(!all(ind)) stop(paste("No parameters for phases", missingPhases, "have been indicated."))
  
  return(phases)
  
}

.checkReplicates = function(replicates, nphases) {
  
  if(any(is.na(replicates))) stop("NAs are not allowed as 'replicates' number.")
  if(is.null(replicates)) return(rep(1, nphases))
  if(length(replicates)==1) return(rep(replicates, nphases)) 
  if(length(replicates)!=nphases) stop("Length of 'replicates' argument must match number of phases.")
  replicates = as.integer(replicates)
  
  if(any(replicates<=0)) stop("All 'replicates' must be positive integers.")
  
  return(replicates)
  
}

.checkBounds =  function(lower, upper, npar) {

  lf = uf = FALSE
  
  if(is.null(lower)) {
    lower = -Inf
    lf = TRUE
  }
  
  if(is.null(upper)) {
    upper = Inf
    uf = TRUE
  }
  
  nl = length(lower)
  nu = length(upper)
  
  if(nl==1) {
    if(npar!=1 & lower!=-Inf) 
      if(!lf) warning("Only one lower bound has been provided, used for all parameters.")
    lower = rep(lower, npar)
    nl = npar
  }
  
  if(nu==1) {
    if(npar!=1 & upper!=Inf) 
      if(!uf) warning("Only one upper bound has been provided, used for all parameters.")    
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

.optPopSize = function(n, selection) floor(0.5*(4 + floor(3*log(n)))/selection)

.checkControl = function(control, method, par, fn, active, skeleton, replicates=1, ...) {
  
  fn = match.fun(fn)

  # specific options must be taken in .optimES.  
  con = list(trace = 0, fnscale = 1, parscale = rep.int(1L, length(which(active))), maxit = NULL, maxgen=NULL,
             abstol = -Inf, reltol = sqrt(.Machine$double.eps), REPORT = 10L, ncores=parallel::detectCores(), 
             alpha=0.05, age.max=1, selection=0.5, step=0.5, nvar=NULL, weights=1, sigma=NULL,
             method=method, aggFn=.weighted.sum, parallel=FALSE, run=NULL, master=NULL, useCV=TRUE,
             convergence=1e-6, stochastic=FALSE, verbose=FALSE, restart.file=NULL)
  
  
  popOpt = .optPopSize(n=length(which(active)), selection=con$selection)
  con$popsize = popOpt
  
  controlDef  = names(con)       # default options
  controlUser = names(control)   # user provided options
  
  con[controlUser] = control
  
  # get unknown options
  unknown = controlUser[!(controlUser %in% controlDef)]
  
  # check for wrong master argument (not a path nor NULL)
  if(!is.character(con$master) & !is.null(con$master)) {
    warning("Control parameter 'master' must be a path or NULL, ignoring user provided value.")
    con$master = NULL
  }
  
  con$ncores = as.integer(con$ncores)
  if(is.null(con$ncores) | is.na(con$ncores)) stop("Control option 'ncores' must be a positive integer.")
  if(is.null(con$ncores)) con$ncores = 1
  if(con$ncores < 1) con$ncores = 1
  
  if(!is.null(con$master) & is.null(con$run)) {
    con$run = tempdir()
    msg = sprintf("Evaluating 'fn' in %s.", con$run)
    message(msg)
  }
  
  if(!is.null(con$run)) {
    if(!dir.exists(con$run)) dir.create(con$run, recursive=TRUE)
  }
  
  # update population size and selection rate
  if(con$popsize < popOpt) {
    warning("'popsize' is too small, using default value.")
    message(sprintf("Population size set to %d (optimal value).", popOpt))
  }
  con$popsize = max(con$popsize, popOpt)
  if(isTRUE(con$parallel)) {
    popOptP = ceiling(con$popsize/con$ncores)*con$ncores
    if(con$popsize != popOptP) {
      message(sprintf("Optimizing 'popsize' to work with %d cores...", con$ncores))
      message(sprintf("   Optimal population size (%d) has been increased to %d.\n", popOpt, popOptP))
    }
      
    con$popsize = popOptP
    con$selection = con$selection*popOpt/popOptP
  }
  
  if(inherits(fn, "objFn")) {
    con$nvar = attr(fn, "nvar")
    con$weights = attr(fn, "weights")
    names(con$weights) = attr(fn, "variables")
    if(con$verbose) {
      
      #print(con$weights) # TO_DO: make nicer output 
    }
  }
  # check number of variables
  xpar = if(missing(skeleton)) par else relist(par, skeleton)
  if(is.null(con$nvar)) con$nvar = length(fn(xpar, ...)) # HERE CHECK
  
  # update maximum number of function evaluations and generations
  if(!is.null(con$maxit) & !is.null(con$maxgen)) 
    warning("'maxit' and 'maxgen' provided, ignoring 'maxit'.")
  if(!is.null(con$maxit) & is.null(con$maxgen)) 
    con$maxgen = ceiling(con$maxit/con$popsize/replicates)
  if(is.null(con$maxit) & is.null(con$maxgen)) con$maxgen = 2000L
  
  con$maxit = con$popsize*con$maxgen

  if(method=="default") {
    
    # check for user-provided variances (sigma)
    if(!is.null(con$sigma) & length(con$sigma)!=length(active)) 
      stop("Vector of variances (sigma) must match parameter length.")
    if(!is.null(con$sigma)) con$sigma = con$sigma[which(active)]
    
    # update and check weights
    if(length(con$weights)==1) con$weights = rep(con$weights, con$nvar)
    if(length(con$weights)!=con$nvar) stop("Vector of weights should match the length of the output of fn.")
    if(any(con$weights<0)) stop("Weights should be positive numbers.")
    if(any(is.na(con$weights))) stop("Weights cannot be NA.")    
    
  }
  
  if(method=="cmaes") con$weights = NULL 
  # aggregation function for global fitness
  con$aggFn = match.fun(con$aggFn)
  
  if(length(unknown)!=0) 
    warning("Unknown control parameters: ", paste0(paste(unknown, collapse = ", ")),".")
  
  return(con)
  
}

.checkConvergence = function(control, nphases) {
  
  maxgen      = control$maxgen
  maxiter     = control$maxiter
  convergence = control$convergence
  
  if(length(maxgen)==1) maxgen = rep(maxgen, nphases)
  if(length(maxiter)==1) maxiter = rep(maxiter, nphases)
  if(length(convergence)==1) convergence = rep(convergence, nphases)
  
  if(!is.null(maxgen) & length(maxgen)!=nphases) 
    stop("'maxgen' length must match number of calibration phases.")
  
  if(!is.null(maxiter) & length(maxiter)!=nphases) 
    stop("'maxiter' length must match number of calibration phases.")
  
  if(!is.null(convergence) & length(convergence)!=nphases) 
    stop("'convergence' length must match number of calibration phases.")
  
  return(list(maxgen=maxgen, maxiter=maxiter, convergence=convergence))
  
}

.setWorkDir = function(run, i) {
  
  if(is.null(run)) return(invisible())
  workDir = file.path(run, paste0("i", i))
  if(!dir.exists(workDir)) dir.create(workDir, recursive=TRUE)
  setwd(workDir)
  return(workDir)
  
}

.getWorkDir = function(run, i) {
  
  if(is.null(run)) return(invisible())
  workDir = file.path(run, paste0("i", i))
  return(workDir)
  
}


.updateControl = function(control, opt, method) {
  # update sigma...
  return(control)
}

.myRowMean = function(x, na.rm=FALSE) {
  if(!inherits(x, "array") & !inherits(x, "matrix")) return(NULL)
  if(length(dim(x)) <  2) return(x)
  if(length(dim(x)) == 2) return(rowMeans(x, na.rm=na.rm))
  if(length(dim(x)) >  2) return(apply(x, MARGIN=seq_along(dim(x)[-1]), FUN=mean, na.rm=na.rm))
}


.mycbind = function(...) {
  
  xall = list(...)
  x = xall[[1]]
  
  if(length(dim(x)) <  2) return(cbind(...))
  
  out = array(dim=c(dim(x), length(xall)))
  
  if(length(dim(x)) == 2) {
    for(i in seq_along(xall)) {
      out[, , i] = xall[[i]]
    }
    return(out)
  }
  
  if(length(dim(x)) == 3) {
    for(i in seq_along(xall)) {
      out[, , , i] = xall[[i]]
    }
    return(out)
  }
  
  if(length(dim(x)) == 4) {
    for(i in seq_along(xall)) {
      out[, , , , i] = xall[[i]]
    }
    return(out)
  }
  
  stop(sprintf("Dimensions higher than %d are not currently supported.", length(dim(x))-1))
  
}


.write_calibrar_dump = function(run, gen, i) {
  if(is.null(run)) return(invisible(NULL))
  dump = file.path(run, "..", "_calibrar.dump")
  if(!dir.exists(dump)) dir.create(dump, recursive = TRUE)
  ipath = file.path(run, sprintf("i%d", i))
  opath = file.path(dump, sprintf("gen%d_i%d", gen, i))
  if(!dir.exists(opath)) dir.create(opath, recursive=TRUE)
  file.copy(from=file.path(ipath, dir(ipath)), to=opath, recursive=TRUE)
  msg = "There was a problem running your function (generation %d, individual %d). Check '%s' for debugging."
  message(sprintf(msg, gen, i, opath))
  return(invisible(NULL))
}

.mylength = function(x, n) {
  if(length(x)!=1) return(x)
  length(x) = n
  x[is.na(x)] = Inf
  return(x)
}


.rbind_fitness = function(x) {
  n = max(sapply(x, FUN=length))
  out = do.call(rbind, lapply(x, FUN=.mylength, n=n))
  return(out)
}



