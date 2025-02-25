#  Author:  John C Nash
#  Date:    Jul 5, 2016 update
# hjn.R -- R implementation of J Nash BASIC HJG.BAS 20160705
hjn = function(par, fn, lower=-Inf, upper=Inf, bdmsk=NULL, control=list(trace=0), ...) {
  
  restart = .restartCalibration(control) # flag: TRUE or FALSE
  
  n = as.integer(length(par))  # number of elements in par vector

  ctrl =  list(trace=0, maxfeval=2000*n, maximize = FALSE, 
               stepsize=1, stepredn=0.1, eps = 1e-07, reltest=100, 
               restart.file=NULL, REPORT=1L, verbose=FALSE)
  
  namc = names(control)
  missnm = paste(namc[!(namc %in% names(ctrl))], collapse=", ")
  if (!all(namc %in% names(ctrl))) 
    stop(sprintf("unknown names in control: %s", missnm))
  ctrl[namc] = control  #
  control = ctrl
  
  if(!is.null(control$maximize) && control$maximize) 
    stop("Do NOT try to maximize with hjn()")
  
  # if (is.null(control$trace)) control$trace = 0 # just in case
  # if (is.null(control$stepsize)) {
  #    stepsize = 1 # initial step size (could put in control())
  # } else { stepsize = control$stepsize }
  # Want stepsize positive or bounds get messed up
  # if (is.null(control$stepredn)) {
  #    stepredn = .1 # defined step reduction (again in control()?)
  # } else { stepredn = control$stepredn }
  # if (is.null(control$maxfeval)) control$maxfeval=2000*n
  # if (is.null(control$eps)) control$eps=1e-07

  # if (is.null(control$reltest)) { reltest = 100 } else {reltest = control$reltest}

  stepsize = control$stepsize
  stepredn = control$stepredn
  steptol = control$eps  
  reltest = control$reltest
  
  # Hooke and Jeeves with bounds and masks
  if (length(upper) == 1) upper = rep(upper, n)
  if (length(lower) == 1) lower = rep(lower, n)
  if (is.null(bdmsk)) { 
      bdmsk = rep(1,n)
      idx = 1:n 
  } else { idx = which(bdmsk != 0) } # define masks
  # if (any(lower >= upper)){
  #     warning("hjn: lower >= upper for some parameters -- set masks")
  #     bdmsk[which(lower >= upper)] = 0
  #     idx = which(bdmsk != 0)
  # }
  if (control$trace > 0) {
    cat("hjn:bdmsk:")
    print(bdmsk)
  }
  nac = length(idx)
  if (any(par < lower) || any(par > upper)) stop("hjn: initial parameters out of bounds")
 
  if(isTRUE(restart)) {
    
    res = .getRestart(control=control)
    opt = res$opt
    trace = res$trace
    
    keepgoing = opt$keepgoing
    
    par       = opt$par
    f         = opt$f
    fmin      = opt$fmin
    fold      = opt$fold 
    pbase     = opt$pbase
    pbest     = opt$pbest
    fcount    = opt$fcount
    ccode     = opt$ccode
    stepsize  = opt$stepsize
    
    msg       = opt$msg
    gen       = opt$gen

    .messageByGen(opt=opt, trace=trace, restart=TRUE, method="hjn")
    
  } else {
    
    keepgoing = TRUE
    pbase = par # base parameter set (fold is function value)
    f = fn(par, ...) # fn at base point
    fmin = fold = f # "best" function so far
    pbest = par # Not really needed 
    fcount = 1 # count function evaluations, compare with maxfeval
    ccode = 1 # start assuming won't get to solution before feval limit
    
    msg = "Not converged."
    
    trace = list()
    gen = 0
    
  }
    
  while (keepgoing) {
    gen = gen + 1
    # exploratory search -- fold stays same for whole set of axes
    if (control$trace > 0) cat("Exploratory move - stepsize = ", stepsize,"\n")
    if (control$trace > 1) {
       cat("p-start:")
       print(par)
    }
    for (jj in 1:nac) { # possibly could do this better in R
       # use unmasked parameters
       j = idx[jj]
       ptmp = par[j]
       doneg = TRUE # assume we will do negative step
       if(ptmp + reltest < upper[j] + reltest) { # Not on upper bound so do pos step 
          par[j] = min(ptmp+stepsize, upper[j])
          if((par[j] + reltest) != (ptmp + reltest)) {
             fcount = fcount + 1
             f = fn(par, ...)
             if(f < fmin) {
                fmin = f
                pbest = par
                doneg = FALSE # only case where we don't do neg
                resetpar = FALSE
             } 
          } 
       } # end not on upper bound
       if(fcount >= control$maxfeval) break
       if(doneg) {
         resetpar = TRUE # going to reset the parameter unless we improve
         if ((ptmp + reltest) > (lower[j] + reltest)) { # can do negative step
            par[j] = max((ptmp - stepsize), lower[j])
            if ((par[j] + reltest) != (ptmp + reltest)) {
               fcount = fcount + 1
               f = fn(par, ...)
               if (f < fmin) {
                  fmin = f
                  pbest = par
                  resetpar = FALSE # don't reset parameter
               } 
            }
         } #  neg step possible
       } # end doneg
       if(resetpar) par[j] = ptmp
    } # end loop on axes
    if (fcount >= control$maxfeval) {
        ccode = 1
        msg = "Function count limit exceeded."
        counts = c('function'=fcount, 'gradient'=NA)
        ans = list(par=pbest, value=fmin, counts=counts, convergence=ccode, message=msg, trace=trace)
        if (control$trace > 0) cat(msg, "\n")
        return(ans)
    }
    if (control$trace > 0) { 
       cat("axial search with stepsize =",stepsize,"  fn value = ",fmin,"  after ",fcount,"  maxfeval =", control$maxfeval,"\n")
    }
    if (fmin < fold) { # success -- do pattern move (control$trace > 0) cat("Pattern move \n")
       if (control$trace > 1) {
          cat("PM from:")
          print(par)
          cat("pbest:")
          print(pbest)
       }
       for (jj in 1:nac) { # Note par is best set of parameters
          j = idx[jj]
          ptmp = 2.0*par[j] - pbase[j]
          if (ptmp > upper[j]) ptmp = upper[j]
          if (ptmp < lower[j]) ptmp = lower[j]
          pbase[j] = par[j]
          par[j] = ptmp 
       }
       fold = fmin
       if (control$trace > 1) {
          cat("PM to:")
          print(par)
       }

    } else { # no success in Axial Search, so reduce stepsize
       if (fcount >= control$maxfeval) {
        ccode = 1
        msg = "Function count limit exceeded."
        counts = c('function'=fcount, 'gradient'=NA)
        ans = list(par=pbest, value=fmin, counts=counts, convergence=ccode, message=msg, trace=trace)
        if (control$trace > 0) cat(msg, "\n")
        return(ans)
       }
       # first check if point changed
       samepoint = identical((par + reltest),(pbase + reltest))
       if (samepoint) { 
          stepsize = stepsize*stepredn
          if (control$trace > 1) cat("Reducing step to ",stepsize,"\n")
          if (stepsize <= steptol) keepgoing = FALSE
          ccode = 0 # successful convergence
       } else { # return to old base point
          if (control$trace > 1) {
             cat("return to base at:")
             print(pbase)
          }
          par = pbase
       }
    }
    if (fcount >= control$maxfeval) {
      ccode = 1
      msg = "Function count limit exceeded."
      counts = c('function'=fcount, 'gradient'=NA)
      ans = list(par=pbest, value=fmin, counts=counts, convergence=ccode, message=msg, trace=trace)
      if (control$trace > 0) cat(msg, "\n")
      return(ans)
    }
    
    opt = list(keepgoing = keepgoing, par=par, f=f, fmin=fmin, fold=fold, 
               pbase=pbase, pbest=pbest, fcount=fcount, ccode=ccode, 
               stepsize=stepsize, msg=msg, gen=gen, control=control)
    
    .createRestartFile(opt=opt, trace=trace, control=control, method="hjn") # (opt, control)
    
    if(control$verbose & opt$gen%%control$REPORT==0) {
      .messageByGen(opt=opt, trace=trace, method="hjn")
    }
    
  } # end keepgoing loop 
  
  if (control$trace > 1) {
    if(identical(pbest, pbase)) {
      cat("pbase = pbest\n") 
    } else { 
      cat("BAD!: pbase != pbest\n") 
    } 
  }
   
  msg = "hjn appears to have converged"
  counts = c('function'=fcount, 'gradient'=NA)
  ans = list(par=pbest, value=fmin, counts=counts, convergence=ccode, message=msg, trace=trace)
  
  return(ans)  
  
}
