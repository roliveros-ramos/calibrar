#  Author:  John C Nash
#  Date:    Dec 6, 2021 update
#
## An R version of the Nash version of Fletcher's Variable
#   Metric minimization -- bounds constrained parameters
# This uses a simple backtracking line search.
#
# Input:
# par  = a vector containing the starting point
# fn = objective function (assumed to be sufficiently
#   differentiable)
# gr = gradient of objective function, provided as a function
#   or the character name of a numerical approximation function
#  lower = vector of lower bounds on parameters
#  upper = vector of upper bounds on parameters
# Note: free parameters outside bounds will be adjusted to
#   bounds unless control$keepinputpar = TRUE.
# bdmsk = control vector for bounds and masks. Parameters
#   for which bdmsk are 1 are unconstrained or 'free', 
#   those with bdmsk 0 are masked i.e., fixed.
# For historical reasons, we use the same array as an
#   indicator that a parameter is at a lower bound (-3) 
#   or upper bound (-1) 
#   # control = list of control parameters
#    maxit = a limit on the number of iterations (default 500)
#    trace = 0 (default) for no output,
#            > 0 for output (bigger => more output)
#    dowarn=TRUE by default. Set FALSE to suppress warnings.
#    eps = a tolerance used for judging small gradient norm
#           (default = 1e-07). See code for usage.
#    maxit = a limit on the gradient evaluations (default
#             500 + 2*n )
#    maxfeval = a limit on the function evaluations (default
#             3000 + 10*n )
#    maximize = TRUE to maximize the function (default FALSE)
#    reltest = 100.0 (default). Additive shift for equality test.
#    stopbadupdate = TRUE (default). Don't stop when steepest
#             descent search point results in failed inverse 
#             Hessian update
#
# Output:
#    A list with components:
#
#   par: The best set of parameters found.
#
#   value: The value of 'fn' corresponding to 'par'.
#
#   counts: A two-element integer vector giving the number of
#     calls to 'fn' and 'gr' respectively. This excludes those calls
#     needed to compute the Hessian, if requested, and any 
#     calls to 'fn' to compute a finite-difference approximation 
#     to the gradient.
#
#   convergence: An integer termination code. 
#      '0' indicates that Rvmmin judges that successful 
#          termination has been obtained.
#       other termination codes are
#          '0' converged, apparently successfully
#          '1' indicates that the maximum iterations 'maxit' or
#              function evaluation count 'maxfeval' was reached.
#          '2' indicates that a point has been found with small
#              gradient norm (< (1 + abs(fmin))*eps*eps )
#          '3' indicates approx. inverse Hessian cannot be updated
#              at steepest descent iteration (i.e., something 
#              very wrong)
#          '20' indicates initial point is infeasible/inadmissible
#          '21' indicates a set of parameters has been tried that
#               are infeasible (function cannot be computed)
#
#   message: A character string giving any additional
#     information returned by the optimizer, or 'NULL'.
#
#   bdmsk: Returned index describing the status of bounds and 
#     masks at the proposed solution. Parameters for which 
#     bdmsk are 1 are unconstrained or 'free', those with 
#     bdmsk 0 are masked i.e., fixed. For historical
#     reasons, we indicate a parameter is at a lower bound
#     using -3 or upper bound using -1.
#
#################################################################
# control defaults
Rvmminb = function(par, fn, gr = NULL, lower = NULL, 
                    upper = NULL, bdmsk = NULL, control = list(), ...) {

  restart = .restartCalibration(control) # flag: TRUE or FALSE
  
  n = as.integer(length(par))  # number of elements in par vector
  maxit = 500 + 2L * n
  maxfeval = 3000 + 10L * n
  ctrl = list(maxit = maxit, maxfeval = maxfeval, maximize = FALSE, 
               trace = 0, eps = 1e-07, dowarn = TRUE, acctol = 0.0001, stepredn=0.2,
               reltest=100.0, stopbadupdate = TRUE)
  namc = names(control)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  ctrl[namc] = control  #
  maxit = ctrl$maxit  #
  maxfeval = ctrl$maxfeval  #
  maximize = ctrl$maximize  # TRUE to maximize the function
  # trace = ctrl$trace  #
  eps = ctrl$eps  #
  acctol = ctrl$acctol # 130125
  dowarn = ctrl$dowarn  #
  stepredn = ctrl$stepredn
  reltest = ctrl$reltest
  stopbadupdate = ctrl$stopbadupdate
  fargs = list(...)  # the ... arguments that are extra function / gradient data
  smallstep = reltest*.Machine$double.xmin # 20230727 fix for neg trystep
  #################################################################
  # check if there are bounds
  if (is.null(lower) || !any(is.finite(lower))) 
    nolower = TRUE
  else nolower = FALSE
  if (is.null(upper) || !any(is.finite(upper))) 
    noupper = TRUE
  else noupper = FALSE
  if (nolower && noupper && all(bdmsk == 1)) { 
    bounds = FALSE
    stop("Do not use Rvmminb() without bounds.")
  } else { bounds = TRUE }
  
  #################################################################
  ## Set working parameters (See CNM Alg 22)
  bvec = par  # copy the parameter vector
  n = length(bvec)  # number of elements in par vector
  ifn = 1  # count function evaluations
  #  stepredn = 0.2  # Step reduction in line search
  #  reltest = 100  # relative equality test
  ceps = .Machine$double.eps * reltest
  dblmax = .Machine$double.xmax  # used to flag bad function
  #############################################
  # gr MUST be provided
  if (is.null(gr)) {  # if gr function is not provided STOP 
    stop("A gradient calculation (analytic or numerical) MUST be provided for Rvmminb")
  }
  mygr = gr
  
  f = try(fn(bvec, ...), silent=FALSE) # Compute the function.
  if (inherits(f,"try-error") | is.na(f) | is.null(f) | is.infinite(f)) {
    msg = "Initial point gives inadmissible function value"
    conv = 20
    ans = list(bvec, dblmax, c(ifn, 0), conv, msg, bdmsk)  #
    names(ans) = c("par", "value", "counts", "convergence", 
                    "message", "bdmsk")
    return(ans)
  }
  g = mygr(bvec, ...)  # Do we need to use try() ?
  
  if(isTRUE(restart)) {
    
    res = .getRestart(control=control)
    opt = res$opt
    trace = res$trace
    
    keepgoing = opt$keepgoing
    bvec      = opt$bvec
    par       = opt$par
    ifn       = opt$ifn
    ig        = opt$ig
    ilast     = opt$ilast
    fmin      = opt$fmin
    f         = opt$f
    g         = opt$g
    oldstep   = opt$oldstep
    conv      = opt$conv
    B         = opt$B
    bdmsk     = opt$bdmsk
    msg       = opt$msg
    
    .messageByGen(opt=opt, trace=trace, restart=TRUE, method="Rvvmin")
    
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
    
    trace = list()
    
  }
  
  gnorm = sqrt(sum(g*g)) ## JN180414 
  if (gnorm < (1 + abs(fmin))*eps*eps ) {
    keepgoing = FALSE
    conv = 2
  }
  
  # START OF MAIN LOOP
  while(keepgoing) { ## main loop -- must remember to break out of it!
    if (ilast == ig) { # reset the approx. inverse hessian B to unit matrix
      B = diag(1, n, n)  # create unit matrix of order n
    }
    # fmin = f # ROR: check
    # par = bvec  # save parameters
    if (!all(is.numeric(g))) {
      g = rep(0, n)  # 110619
      cat("zeroing gradient because of failure\n")
    }
    c = g  # save gradient
    ## Bounds and masks adjustment of gradient ##
    ## current version with looping -- later try to vectorize
    for (i in 1:n) {
      if ((bdmsk[i] == 0)) {
        g[i] = 0
      }
      else {
        if (bdmsk[i] == 1) {
        }
        else {
          if ((bdmsk[i] + 2) * g[i] < 0) {
            g[i] = 0  # active mask or constraint
          }
          else {
            bdmsk[i] = 1  # freeing parameter i
          }
        }
      }
    }  # end masking loop on i
    ## end bounds and masks adjustment of gradient
    t = as.vector(-B %*% g)  # compute search direction
    if (!all(is.numeric(t))) 
      t = rep(0, n)  # 110619
    t[which(bdmsk <= 0)] = 0  # apply masks and box constraints
    gradproj = sum(t * g)  # gradient projection
    accpoint = FALSE  # Need this BEFORE gradproj test
    
    # ROR: here we compute trystep only once
    trystep = ifelse(t<0, (lower - par)/t, (upper - par)/t)
    trystep[(bdmsk != 1) | (t == 0)] = NA
    
    if (is.nan(gradproj)) {
      warning("gradproj Nan")
      gradproj = 0  # force null
    }
    if (gradproj <= 0) {
      # Must be going downhill OR be converged
      ########################################################
      ####      Backtrack only Line search                ####
      changed = TRUE  # Need to set so loop will start
      steplength = oldstep # 131202 - 1 seems best value (Newton step)
      
      while (changed && (!accpoint)) {
        # We seek a lower point, but must change parameters too
        # Box constraint -- adjust step length for free parameters
        # ROR: this is the only update needed, if so (trystep is fix, steplength is always lower)
        steplength = min(steplength, trystep, na.rm=TRUE)
        if(steplength < smallstep) steplength = 0 # force break for neg step
        
        # for (i in 1:n) { # loop on parameters -- vectorize?
        #   if ((bdmsk[i] == 1) && (t[i] != 0)) {
        #     # only concerned with free parameters and non-zero search dimension
        #     if (t[i] < 0) {
        #       # going down. Look at lower bound
        #       trystep = (lower[i] - par[i])/t[i]  # t[i] < 0 so this is positive
        #     } else {
        #       # going up, check upper bound
        #       trystep = (upper[i] - par[i])/t[i]  # t[i] > 0 so this is positive
        #     }
        #     steplength = min(steplength, trystep)  # reduce as necessary
        #     if (steplength < smallstep) steplength = 0 # force break for neg step
        #   }  # end steplength reduction
        # }  # end loop on i to reduce step length
        # end box constraint adjustment of step length
        
        bvec = par + steplength * t
        changed = (!identical((bvec + reltest), (par + reltest)) )
        if (changed) {
          # compute new step, if possible
          f = try(fn(bvec, ...))
          if (inherits(f, "try-error")) f = .Machine$double.xmax
          if (maximize) f = -f
          ifn = ifn + 1
          if (ifn > maxfeval) {
            msg = "Too many function evaluations"
            if (dowarn) warning(msg)
            conv = 1
            changed = FALSE
            keepgoing = FALSE
            break # without saving parameters
          }
          if (is.infinite(f)) f = .Machine$double.xmax
          if (is.na(f) | is.null(f) ) {
            msg='Function is not calculable at an intermediate point'
            #  stop('f is NA')
            conv = 21
            f = dblmax  # try big function to escape
            keepgoing = FALSE
            break
          }
          accpoint = (f < fmin + gradproj * steplength * acctol) # NOTE: < not <=
          if (! accpoint) {
            steplength = steplength * stepredn
          }
        } # end changed
        else { # NOT changed in step reduction
        }
      }  # end while ((f >= fmin) && changed )
    }  # end if gradproj<0
    
    if (accpoint) {
      fmin = f # remember to save the value 150112
      # matrix update if acceptable point.
      for (i in 1:n) { ## Reactivate constraints?
        if (bdmsk[i] == 1) { # only interested in free parameters
          # make sure < not <= below to avoid Inf comparisons
          if ((bvec[i] - lower[i]) < ceps * (abs(lower[i]) + 1)) {
            # are we near or lower than lower bd
            bdmsk[i] = -3
          }  # end lower bd reactivate
          if ((upper[i] - bvec[i]) < ceps * (abs(upper[i]) + 1)) {
            # are we near or above upper bd
            bdmsk[i] = -1
          }  # end lower bd reactivate
        }  # end test on free params
      }  # end reactivate constraints loop
      test = try(g = mygr(bvec, ...), silent = FALSE)
      if (inherits(test, "try-error")) stop("Bad gradient!!")
      if (any(is.nan(g))) stop("NaN in gradient")
      ig = ig + 1
      if (maximize) g = -g
      if (ig > maxit) {
        keepgoing = FALSE
        msg = "Too many gradient evaluations"
        if (dowarn) warning(msg)
        conv = 1
        break
      }
      par = bvec # save parameters since point acceptable
      ## ERROR!!        g[which(bdmsk <= 0)] = 0  # adjust for active mask or constraint 
      if (bounds) 
      { ## Bounds and masks adjustment of gradient ##
        ## first try with looping -- later try to vectorize
        for (i in 1:n) {
          if ((bdmsk[i] == 0)) {
            # masked, so gradient component is zero
            g[i] = 0
          }
          else {
            if (bdmsk[i] == 1) {
            }
            else {
              if ((bdmsk[i] + 2) * g[i] < 0) {
                # test for -ve gradient at upper bound, +ve at lower bound
                g[i] = 0  # active mask or constraint and zero gradient component
              }
              else {
                bdmsk[i] = 1  # freeing parameter i
              }
            }
          }
        }  # end masking loop on i
      }  # end if bounds
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
        # May be able to be more efficient below -- need to use
        #   outer function
        B = B - (outer(t, y) + outer(y, t) - D2 * outer(t, t))/D1
      }
      else {
        if (ig == ilast+1) {
          if (stopbadupdate && ! accpoint) keepgoing=FALSE # stop on update failure for s.d. search
          conv = 3
        }
        ilast = ig  # note gradient evaluation when update failed
      }  # D1 > 0 test
    } # end if accpoint
    else { # no acceptable point
      if ( (ig == ilast) || (abs(gradproj) < (1 + abs(fmin))*eps*eps ) ) { # remove ig > 2
        # we reset to gradient and did new linesearch
        keepgoing = FALSE  # no progress possible
        if (conv < 0) { # conv == -1 is used to indicate it is not set
          conv = 0
        }
        msg = "Converged"
      } # end ig == ilast
      else {
        ilast = ig  # reset to gradient search
      }  # end else ig != ilast
    }  # end else no accpoint

    opt = list(keepgoing=keepgoing, bvec=bvec, par=par, ifn=ifn, ig=ig, 
               ilast=ilast, fmin=fmin, f=f, g=g, oldstep=oldstep,
               conv=conv, B=B, bdmsk=bdmsk, msg=msg)
              
    .createRestartFile(opt=opt, trace=trace, control=control, method="Rvvmin") # (opt, control)
        
  }  # end main loop  (while keepgoing)
  
  msg = "Rvmminb appears to have converged"
  counts = c('function'=ifn, 'gradient'=ig)
  ans = list(par=par, value=fmin, counts=counts, convergence=conv, message=msg, trace=trace)
  
  return(ans) 
  
}  ## end of Rvmminb
