
# New Rvmmin --------------------------------------------------------------

.Rvmminb2 = function(par, fn, gr = NULL, lower = NULL, upper = NULL, control = list(), bdmsk=NULL, ...) {
  # control defaults
  
  restart = .restartCalibration(control) # flag: TRUE or FALSE
  
  n = as.integer(length(par))  # number of elements in par vector
  bdmsk = rep(1L, n) # all parameters are active, it's controlled outside
  ctrl = list(maxit = 500 + 2L * n, maxfeval = 3000 + 10L * n, trace = 0, eps = 1e-07, 
              dowarn = TRUE, acctol = 0.0001, stepredn=0.2,
              reltest=100, stopbadupdate = TRUE)
  namc = names(control)
  
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  
  ctrl[namc] = control  #
  maxit = ctrl$maxit  #
  maxfeval = ctrl$maxfeval  #
  trace = ctrl$trace  #
  eps = ctrl$eps  #
  acctol = ctrl$acctol # 130125
  dowarn = ctrl$dowarn  #
  stepredn = ctrl$stepredn
  reltest = ctrl$reltest
  stopbadupdate = ctrl$stopbadupdate
  fargs = list(...)  # the ... arguments that are extra function / gradient data
  # check if there are bounds
  nolower = if(is.null(lower) || !any(is.finite(lower))) TRUE else FALSE
  noupper = if(is.null(upper) || !any(is.finite(upper))) TRUE else FALSE
  bounds = if(nolower && noupper) FALSE else TRUE 
  if(!bounds) stop("Do not use Rvmminb() without bounds.")
  
  ## Set working parameters (See CNM Alg 22)
  ceps = .Machine$double.eps * reltest
  dblmax = .Machine$double.xmax  # used to flag bad function
  
  # Assume bounds already checked 150108
  f = try(fn(par, ...), silent=TRUE) # Compute the function.
  
  if((inherits(f, "try-error")) | is.na(f) | is.null(f) | is.infinite(f)) {
    msg = "Initial point gives inadmissible function value"
    conv = 20
    ans = list(par=par, value=dblmax, counts=c(1, 0), convergence=conv, message=msg)
    return(ans)
  }
  
  # gr MUST be provided
  if(is.null(gr)) {  # if gr function is not provided STOP 
    stop("A gradient calculation (analytic or numerical) MUST be provided for Rvmminb")
  }
  mygr = gr
  
  if(isTRUE(restart)) {
    
    res = .getRestart(control=control)
    
    keepgoing = res$keepgoing
    bvec      = res$bvec
    par       = res$par
    ifn       = res$ifn
    ig        = res$ig
    ilast     = res$ilast
    fmin      = res$fmin
    f         = res$f
    g         = res$g
    oldstep   = res$oldstep
    conv      = res$conv
    B         = res$B
    bdmsk     = res$bdmsk
    msg       = res$msg
    # .messageByGen(opt, trace, restart=TRUE)
    
  } else {
    
    keepgoing = TRUE  # to ensure loop continues until we are finished
    # initialization
    bvec = par  # copy the parameter vector
    ifn = 1  # count function evaluations # RESTART
    ig = 1  # count gradient evaluations # RESTART
    ilast = ig  # last time we used gradient as search direction # RESTART
    fmin = f  # needed for numerical gradients # RESTART
    g = mygr(bvec, ...)  # Do we need to use try() ?
    oldstep = 1
    conv = -1
    fmin = f # Initialize fmin
  }
  
  gnorm = sqrt(sum(g*g)) ## JN180414
  if (gnorm < (1 + abs(fmin))*eps*eps ) {
    keepgoing = FALSE
    conv = 2
  }
  
  ### TRACE, to update after restart
  BDMSK = matrix(nrow=maxit, ncol=n)
  jj = 1
  
  g[which(bdmsk == 0)] = 0 # fixed points, not needed now.
  
  ## START OF THE OPTIMIZATION
  while(keepgoing) { ## main loop -- must remember to break out of it!
    
    if(ilast == ig) { # reset the approx. inverse hessian B to unit matrix
      B = diag(1, n, n)  # create unit matrix of order n
    }
    
    # fmin = f # ROR: This is an error, f may come from a non acceptable point.
    # par = bvec  # save parameters # ROR: This is an error too.
    if (!all(is.numeric(g))) {
      g = rep(0, n)  # 110619
      cat("zeroing gradient because of failure\n")
    }
    c = g  # save gradient
    ## Bounds and masks adjustment of gradient ##
    t = as.vector(-B %*% g)  # compute search direction
    if (!all(is.numeric(t))) t = rep(0, n)  # 110619
    t[which(bdmsk <= 0)] = 0  # apply masks and box constraints
    gradproj = sum(t * g)  # gradient projection
    
    # ROR: here we compute trystep only once
    trystep = ifelse(t<0, (lower - par)/t, (upper - par)/t)
    trystep[(bdmsk != 1) | (t == 0)] = NA
    
    accpoint = FALSE  # Need this BEFORE gradproj test
    if(is.nan(gradproj)) {
      warning("gradproj NaN")
      gradproj = 0  # force null
    }
    if (gradproj <= 0) {
      # Must be going downhill OR be converged
      ####      Backtrack only Line search                ####
      changed = TRUE  # Need to set so loop will start
      steplength = oldstep # 131202 - 1 seems best value (Newton step)
      
      while ((f >= fmin) && changed && (!accpoint)) {
        # We seek a lower point, but must change parameters too
        # Box constraint -- adjust step length for free parameters
        # ROR: this is the only update needed, if so (trystep is fix, steplength is always lower)
        steplength = min(steplength, trystep, na.rm=TRUE)
        
        bvec = par + steplength * t
        changed = !identical((bvec + reltest), (par + reltest))
        
        if(changed) {
          # compute new step, if possible
          f = try(fn(bvec, ...))
          if(inherits(f, "try-error")) f = .Machine$double.xmax
          ifn = ifn + 1
          if (ifn > maxfeval) {
            msg = "Too many function evaluations"
            if(dowarn) warning(msg)
            conv = 1
            changed = FALSE
            keepgoing = FALSE
            break # without saving parameters
          }
          if(is.infinite(f)) f = .Machine$double.xmax
          if(is.na(f) | is.null(f) ) {
            msg='Function is not calculable at an intermediate point'
            conv = 21
            f = dblmax  # try big function to escape
            keepgoing = FALSE
            break
          }
          if (f < fmin) {
            # We have a lower point. Is it 'low enough' i.e. acceptable
            accpoint = (f <= fmin + gradproj * steplength * acctol)
          } else {
            steplength = steplength * stepredn
          }
        } # end changed
        
      }  # end while ((f >= fmin) && changed )
      
    }  # end if gradproj<0
    
    if(accpoint) {
      fmin = f # remember to save the value 150112
      # matrix update if acceptable point.
      
      # for(i in 1:n) { ## Reactivate constraints?
      #   if (bdmsk[i] == 1) { # only interested in free parameters
      # 
      #   }  # end test on free params
      # }  # end reactivate constraints loop
    
      g = try(mygr(bvec, ...), silent = TRUE) 
      if (inherits(g, "try-error")) stop("Bad gradient!!")
      if(any(is.nan(g))) stop("NaN in gradient")
      ig = ig + 1
      if(ig > maxit) {
        keepgoing = FALSE
        msg = "Too many gradient evaluations"
        if(dowarn) warning(msg)
        conv = 1
        break
      }
      par = bvec # save parameters since point acceptable
      
      ## current version with looping -- later try to vectorize
      for (i in 1:n) {
        if((bdmsk[i] == 0)) {
          g[i] = 0
        } else {
          if(bdmsk[i] == 1) { # currently active
            # make sure < not <= below to avoid Inf comparisons
            # ROR: compare to range to make it invariant up to a constant
            if ((bvec[i] - lower[i]) < ceps * (abs(lower[i]) + 1)) {
              # are we near or lower than lower bd
              bdmsk[i] = -3
            }  # end lower bd reactivate
            if ((upper[i] - bvec[i]) < ceps * (abs(upper[i]) + 1)) {
              # are we near or above upper bd
              bdmsk[i] = -1
            }  # end lower bd reactivate
          } else { # currently inactive, to reactivate?
            if((bdmsk[i] + 2) * g[i] < 0) {
              # still going in the bad direction, stop it!
              # g[i] = 0  # active mask or constraint
            } else {
              # going in the right direction, activate it!
              bdmsk[i] = 1  # freeing parameter i
            }
          }
        }
      }  # end masking loop on i
      ## end bounds and masks adjustment of gradient
      
      # ROR: CHECK! should we remove from the gnorm values in the border?
      g[which(bdmsk <= 0)] = 0  # adjust for active mask or constraint 
      # ind = which(bdmsk <= 0) 
      gnorm = sqrt(sum(g*g)) ## JN131202 
      if (gnorm < (1 + abs(fmin))*eps*eps ) {
        keepgoing = FALSE
        conv = 2
        break
      }
      ## 150107 check on breakout
      ## if (! keepgoing) stop("break with small gnorm failed")
      t = as.vector(steplength * t)
      c = as.vector(g - c)
      D1 = sum(t * c)
      if (D1 > 0) {
        y = as.vector(crossprod(B, c))
        D2 = as.double(1+crossprod(c,y)/D1)  
        # as.double because D2 is a 1 by 1 matrix otherwise
        B = B - (outer(t, y) + outer(y, t) - D2 * outer(t, t))/D1
      } else {
        if (ig == ilast+1) {
          if(stopbadupdate && !accpoint) keepgoing=FALSE # stop on update failure for s.d. search
          conv = 3
        }
        ilast = ig  # note gradient evaluation when update failed
      }  # D1 > 0 test
      
    } else { # no acceptable point
      
      if( (ig == ilast) || (abs(gradproj) < (1 + abs(fmin))*eps*eps) ) { # remove ig > 2
        # we reset to gradient and did new linesearch
        keepgoing = FALSE  # no progress possible
        if(conv < 0) conv = 0 # conv == -1 is used to indicate it is not set
        msg = "Converged."
      } else {
        ilast = ig  # reset to gradient search
      }  # end else ig != ilast
      
    }  # end else no accpoint
    
    .createRestartFile_Rvmmin() # (opt, control)
    
    BDMSK[jj, ] = bdmsk
    jj = jj + 1
    
  }  # end main loop  (while keepgoing)
  # END OF OPTIMIZATION
  
  BDMSK = BDMSK[1:(jj-1), ]
  
  msg = "Rvmminb appears to have converged"
  counts = c(ifn, ig)
  names(counts) = c("function", "gradient")
  ans = list(par=par, value=fmin, counts=counts, convergence=conv, message=msg, trace=list(bdmsk=BDMSK))
  
  return(ans)
  
}  ## end of Rvmminb
