# less verbose version of Rvmmin::Rvmmin
Rvmmin2 = function(par, fn, gr = NULL, lower = NULL, 
                   upper = NULL, bdmsk = NULL, control = list(), ...) {
  #
  #  Author:  John C Nash
  #  Date:  April 2, 2009; revised July 28, 2009, May 21, 2012, Jan 8, 2015
  #
  ## An R version of the Nash version of Fletcher's Variable
  #   Metric minimization
  # This uses a simple backtracking line search.
  # This is a driver for both constrained and unconstrained routines.
  #
  # Input:
  # par  = a vector containing the starting point
  # fn = objective function (assumed to be sufficiently
  #   differentiable)
  # gr = gradient of objective function
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
  # control = list of control parameters
  #    maxit = a limit on the number of iterations (default 500)
  #    maximize = TRUE to maximize the function (default FALSE)
  #    trace = 0 (default) for no output,
  #            > 0 for output (bigger => more output)
  #    dowarn=TRUE by default. Set FALSE to suppress warnings.
  #    checkgrad = TRUE by default. Check analytic gradient 
  #            against numDeriv results.
  #    checkbounds = TRUE by default. Check parameters and bounds
  #            for addmissible bounds and feasible start.
  #    keepinputpar = FALSE if Rvmmin will move out-of-bound parameters
  #            to nearest bound WITH WARNING. This is default behaviour.
  #            If TRUE, stop execution. That is, do not allow starts to
  #            be modified by bmchk(). This is only in Rvmmin(), not
  #            Rvmminu() or Rvmminb().
  #    eps = a tolerance used for judging small gradient norm
  #           (default = 1e-07). See code for usage.
  #    stepredn = 0.2 (default). Step reduction factor for backtrack
  #             line search
  #    reltest = 100.0 (default). Additive shift for equality test.
  #    stopbadupdate = FALSE (default). Don't stop when steepest
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
  #   convergence: An integer code. '0' indicates that Rvmmin judges
  #     that successful termination has been obtained.
  #          convergence codes are
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
  npar <- length(par) # number of parameters
  # control defaults -- idea from spg
  if (is.null(control$trace)) control$trace=0
  # check if there are bounds
  if (is.null(lower) || !any(is.finite(lower))) 
    nolower = TRUE # no lower bounds
  else nolower = FALSE
  if (is.null(upper) || !any(is.finite(upper))) 
    noupper = TRUE # no upper bounds
  else noupper = FALSE
  if (nolower && noupper && all(bdmsk == 1)) 
    bounds = FALSE
  else bounds = TRUE
  if (control$trace > 1) 
    cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, 
        " bounds = ", bounds, "\n")
  if (is.null(gr)) {
    gr <- "grfwd" # use forward gradient approximation if no gradient code provided
    if (control$trace > 0) cat("WARNING: forward gradient approximation being used\n")
  } else {
    if (is.character(gr)) { # assume numerical gradient
      if (control$trace > 0) cat("WARNING: using gradient approximation '",gr,"'\n")
    } else { # analytic gradient, so check if requested
      if (is.null(control$checkgrad)) control$checkgrad <- TRUE
      if (control$checkgrad) { # check gradient
        testgrad <- grchk(par, fn, gr, trace=control$trace, ...)
        if (! testgrad) warning("Gradient code for Rvmmin may be faulty - check it!")
      }
    } # end else
  }
  control$checkgrad<-NULL # to avoid problems in subsidiary routines
  if (is.null(control$dowarn)) control$dowarn<-TRUE
  #############################################
  if (bounds) { 
    if (is.null(control$checkbounds)) { control$checkbounds <- TRUE }
    if (is.null(bdmsk)) { bdmsk <- rep(1, npar) } # ensure we have bdmsk
    if ((length(lower) == 1) && (npar > 1) ) lower <- rep(lower, npar) 
    if ((length(upper) == 1) && (npar > 1) ) upper <- rep(upper, npar) # fix 150604
    if (any(is.infinite(lower))) lower[which(is.infinite(lower))] <- -.Machine$double.xmax
    if (any(is.infinite(upper))) upper[which(is.infinite(upper))] <-  .Machine$double.xmax
    ### Check bounds feasible
    if (control$checkbounds) {
      btest <- bmchk(par, lower = lower, upper = upper, bdmsk = bdmsk, 
                     trace = control$trace)
      if (!btest$admissible) 
        stop("Inadmissible bounds: one or more lower > upper")
      if (btest$parchanged) {
        if (is.null(control$keepinputpar) || ! control$keepinputpar) { 
          warning("Parameter out of bounds has been moved to nearest bound")
          control$keepinputpar <- NULL # avoid problems in subsidiary routines
          par <- btest$bvec # save the changed parameters             
        } else stop("Parameter out of bounds")
      }
    }
    nolower <- btest$nolower
    noupper <- btest$noupper
    bounds <- btest$bounds
    bdmsk <- btest$bdmsk  # change bdmsk to values set in bmchk
    if (control$trace > 3) {
      cat("Adjusted bdmsk vector:")
      print(bdmsk)
    }
    lower <- btest$lower
    upper <- btest$upper
    control$checkbounds<-NULL # to avoid problems in subsidiary routines
    ############## end bounds check #############
    ans <- .Rvmminb(par, fn, gr, lower = lower, 
                   upper = upper, bdmsk = bdmsk, control = control, ...)
  } else {
    ans <- Rvmminu(par, fn, gr, control = control, ...)
  } #   return(ans) 
}  ## end of Rvmmin


.Rvmminb = function(par, fn, gr = NULL, lower = NULL, 
                     upper = NULL, bdmsk = NULL, control = list(), ...) {
  #
  #  Author:  John C Nash
  #  Date:  , April 14, 2018 -- update to fix infinite loop
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
  n <- as.integer(length(par))  # number of elements in par vector
  maxit <- 500 + 2L * n
  maxfeval <- 3000 + 10L * n
  ctrl <- list(maxit = maxit, maxfeval = maxfeval, maximize = FALSE, 
               trace = 0, eps = 1e-07, dowarn = TRUE, acctol = 0.0001, stepredn=0.2,
               reltest=100.0, stopbadupdate = TRUE)
  namc <- names(control)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  ctrl[namc] <- control  #
  maxit <- ctrl$maxit  #
  maxfeval <- ctrl$maxfeval  #
  maximize <- ctrl$maximize  # TRUE to maximize the function
  trace <- ctrl$trace  #
  eps <- ctrl$eps  #
  acctol <- ctrl$acctol # 130125
  dowarn <- ctrl$dowarn  #
  stepredn <- ctrl$stepredn
  reltest <- ctrl$reltest
  stopbadupdate <- ctrl$stopbadupdate
  fargs <- list(...)  # the ... arguments that are extra function / gradient data
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
  
  # cat("trace=",trace,"\n") # This is not picked up -- Why? It does.
  
  if (trace > 1) {
    cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, 
        " bounds = ", bounds, "\n") 
  }
  #################################################################
  ## Set working parameters (See CNM Alg 22)
  if (trace > 0) {
    cat("Rvmminb -- J C Nash 2009-2015 - an R implementation of Alg 21\n") }
  bvec <- par  # copy the parameter vector
  n <- length(bvec)  # number of elements in par vector
  if (trace > 0) {
    cat("Problem of size n=", n, "  Dot arguments:\n")
    print(fargs)
  }
  ifn <- 1  # count function evaluations
  stepredn <- 0.2  # Step reduction in line search
  reltest <- 100  # relative equality test
  ceps <- .Machine$double.eps * reltest
  dblmax <- .Machine$double.xmax  # used to flag bad function
  #############################################
  # gr MUST be provided
  if (is.null(gr)) {  # if gr function is not provided STOP 
    stop("A gradient calculation (analytic or numerical) MUST be provided for Rvmminb")
  }
  if (is.character(gr)) { # assume numerical gradient
    # Convert string to function call, assuming it is a numerical gradient function
    if (trace > 0) cat("WARNING: using gradient approximation '",gr,"'\n")
    mygr<-function(par=par, userfn=fn, ...){
      do.call(gr, list(par, userfn, ...))
    }
  } else { 
    mygr<-gr 
  } # end else
  ############# end test gr ####################
  # Assume bounds already checked 150108
  f<-try(fn(bvec, ...), silent=TRUE) # Compute the function.
  if ((class(f) == "try-error") | is.na(f) | is.null(f) | is.infinite(f)) {
    msg <- "Initial point gives inadmissible function value"
    conv <- 20
    if (trace > 0) 
      cat(msg, "\n") # change NA to dblmax 110524
    ans <- list(bvec, dblmax, c(ifn, 0), conv, msg, bdmsk)  #
    names(ans) <- c("par", "value", "counts", "convergence", 
                    "message", "bdmsk")
    return(ans)
  }
  if (maximize) f <- -f
  if (trace > 0) cat("Initial fn=", f, "\n")
  if (trace > 2) print(bvec)
  keepgoing <- TRUE  # to ensure loop continues until we are finished
  ig <- 1  # count gradient evaluations
  ilast <- ig  # last time we used gradient as search direction
  fmin <- f  # needed for numerical gradients
  g <- mygr(bvec, ...)  # Do we need to use try() ?
  if (maximize) g <- -g
  if (trace > 2) {
    cat("g:")
    print(g)
  }
  oldstep <- 1
  conv <- -1
  gnorm <- sqrt(sum(g*g)) ## JN180414 
  if (trace > 0) cat("ig=",ig,"  gnorm=",gnorm,"  ")
  if (gnorm < (1 + abs(fmin))*eps*eps ) {
    if (trace > 1) cat("Small gradient norm\n")
    keepgoing <- FALSE
    conv <- 2
  }
  while (keepgoing) { ## main loop -- must remember to break out of it!
    if (ilast == ig) { # reset the approx. inverse hessian B to unit matrix
      B <- diag(1, n, n)  # create unit matrix of order n
      if (trace > 1) cat("Reset Inv. Hessian approx at ilast = ", ilast, "\n")
    }
    fmin <- f
    if (trace > 0) cat(" ", ifn, " ", ig, " ", fmin, "\n")
    par <- bvec  # save parameters
    if (!all(is.numeric(g))) {
      g <- rep(0, n)  # 110619
      cat("zeroing gradient because of failure\n")
    }
    c <- g  # save gradient
    ## Bounds and masks adjustment of gradient ##
    ## current version with looping -- later try to vectorize
    ##         if (bounds) 
    if (trace > 3) {
      cat("bdmsk:")
      print(bdmsk)
    }
    for (i in 1:n) {
      if ((bdmsk[i] == 0)) {
        g[i] <- 0
      }
      else {
        if (bdmsk[i] == 1) {
          if (trace > 2) 
            cat("Parameter ", i, " is free\n")
        }
        else {
          if ((bdmsk[i] + 2) * g[i] < 0) {
            g[i] <- 0  # active mask or constraint
          }
          else {
            bdmsk[i] <- 1  # freeing parameter i
            if (trace > 1) 
              cat("freeing parameter ", i, "\n")
          }
        }
      }
    }  # end masking loop on i
    if (trace > 3) {
      cat("bdmsk adj:")
      print(bdmsk)
      cat("proj-g:")
      print(g)
    }
    ## end bounds and masks adjustment of gradient
    ###    }  # if bounds
    t <- as.vector(-B %*% g)  # compute search direction
    if (!all(is.numeric(t))) 
      t <- rep(0, n)  # 110619
    if (trace > 2) {
      cat("t:")
      print(t)
    }
    t[which(bdmsk <= 0)] <- 0  # apply masks and box constraints
    if (trace > 2) {
      cat("adj-t:")
      print(t)
    }
    gradproj <- sum(t * g)  # gradient projection
    if (trace > 1) 
      cat("Gradproj =", gradproj, "\n")
    accpoint <- FALSE  # Need this BEFORE gradproj test
    if (is.nan(gradproj)) {
      warning("gradproj Nan")
      gradproj <- 0  # force null
    }
    if (gradproj <= 0) {
      # Must be going downhill OR be converged
      ########################################################
      ####      Backtrack only Line search                ####
      changed <- TRUE  # Need to set so loop will start
      steplength <- oldstep # 131202 - 1 seems best value (Newton step)
      while ((f >= fmin) && changed && (!accpoint)) {
        # We seek a lower point, but must change parameters too
        ###if (bounds) { # MUST have bounds in Rvmminb
        # Box constraint -- adjust step length for free parameters
        for (i in 1:n) { # loop on parameters -- vectorize??
          if ((bdmsk[i] == 1) && (t[i] != 0)) {
            # only concerned with free parameters and non-zero search dimension
            if (t[i] < 0) {
              # going down. Look at lower bound
              trystep <- (lower[i] - par[i])/t[i]  # t[i] < 0 so this is positive
            } else {
              # going up, check upper bound
              trystep <- (upper[i] - par[i])/t[i]  # t[i] > 0 so this is positive
            }
            if (trace > 2) cat("steplength, trystep:", steplength, trystep, "\n")
            steplength <- min(steplength, trystep)  # reduce as necessary
          }  # end steplength reduction
        }  # end loop on i to reduce step length
        # end box constraint adjustment of step length
        if (trace > 1) cat("reset steplength=", steplength, "\n")
        ###  }  # end if bounds
        # end box constraint adjustment of step length
        bvec <- par + steplength * t
        if (trace > 2) {
          cat("new bvec:")
          print(bvec)
        }
        changed <- (!identical((bvec + reltest), (par + reltest)) )
        if (trace > 2) cat("changed =",changed,"\n")
        if (changed) {
          # compute new step, if possible
          f <- try(fn(bvec, ...))
          if (class(f) == "try-error") f <- .Machine$double.xmax
          if (maximize) f <- -f
          if (trace > 2) cat("New f=",f," lower = ",(f < fmin),"\n")
          ifn <- ifn + 1
          if (ifn > maxfeval) {
            msg <- "Too many function evaluations"
            if (dowarn) warning(msg)
            conv <- 1
            changed <- FALSE
            keepgoing <- FALSE
            break # without saving parameters
          }
          if (is.infinite(f)) f <- .Machine$double.xmax
          if (is.na(f) | is.null(f) ) {
            if (trace > 2) {
              cat("Function is not calculable at intermediate bvec:")
              print(bvec)
            }
            msg='Function is not calculable at an intermediate point'
            #  stop('f is NA')
            conv <- 21
            f <- dblmax  # try big function to escape
            keepgoing <- FALSE
            break
          }
          if (f < fmin) {
            # We have a lower point. Is it 'low enough' i.e.,
            #   acceptable
            accpoint <- (f <= fmin + gradproj * steplength * acctol)
            if (trace > 2) cat("accpoint = ", accpoint,"\n")
          }
          else {
            steplength <- steplength * stepredn
            if (trace > 0) cat("*")
          }
        } # end changed
        else { # NOT changed in step reduction
          if (trace > 1) cat("Unchanged in step redn \n")
        }
      }  # end while ((f >= fmin) && changed )
    }  # end if gradproj<0
    if (accpoint) {
      fmin <- f # remember to save the value 150112
      # matrix update if acceptable point.
      ### if (bounds) {
      for (i in 1:n) { ## Reactivate constraints?
        if (bdmsk[i] == 1) { # only interested in free parameters
          # make sure < not <= below to avoid Inf comparisons
          if ((bvec[i] - lower[i]) < ceps * (abs(lower[i]) + 1)) {
            # are we near or lower than lower bd
            if (trace > 2) cat("(re)activate lower bd ", i, " at ", lower[i], "\n")
            bdmsk[i] <- -3
          }  # end lower bd reactivate
          if ((upper[i] - bvec[i]) < ceps * (abs(upper[i]) + 1)) {
            # are we near or above upper bd
            if (trace > 2) cat("(re)activate upper bd ", i," at ", upper[i], "\n")
            bdmsk[i] <- -1
          }  # end lower bd reactivate
        }  # end test on free params
      }  # end reactivate constraints loop
      ###   }  # if bounds
      test <- try(g <- mygr(bvec, ...), silent = TRUE) 
      if (class(test) == "try-error") stop("Bad gradient!!")
      if (any(is.nan(g))) stop("NaN in gradient")
      ig <- ig + 1
      if (maximize) g <- -g
      if (ig > maxit) {
        keepgoing = FALSE
        msg = "Too many gradient evaluations"
        if (dowarn) warning(msg)
        conv <- 1
        break
      }
      par <- bvec # save parameters since point acceptable
      g[which(bdmsk <= 0)] <- 0  # adjust for active mask or constraint 
      gnorm <- sqrt(sum(g*g)) ## JN131202 
      if (trace > 0) cat("ig=",ig,"  gnorm=",gnorm,"  ")
      if (gnorm < (1 + abs(fmin))*eps*eps ) {
        if (trace > 1) cat("Small gradient norm\n")
        keepgoing <- FALSE
        conv <- 2
        break
      }
      ## 150107 check on breakout
      ## if (! keepgoing) stop("break with small gnorm failed")
      t <- as.vector(steplength * t)
      c <- as.vector(g - c)
      D1 <- sum(t * c)
      if (D1 > 0) {
        y <- as.vector(crossprod(B, c))
        D2 <- as.double(1+crossprod(c,y)/D1)  
        # as.double because D2 is a 1 by 1 matrix otherwise
        # May be able to be more efficient below -- need to use
        #   outer function
        B <- B - (outer(t, y) + outer(y, t) - D2 * outer(t, t))/D1
      }
      else {
        if (trace > 0) 
          cat("UPDATE NOT POSSIBLE: ilast, ig",ilast, ig,"\n")
        if (ig == ilast+1) {
          if (stopbadupdate && ! accpoint) keepgoing=FALSE # stop on update failure for s.d. search
          if (trace > 2) cat("keepgoing = ",keepgoing,"\n")
          conv <- 3
        }
        ilast <- ig  # note gradient evaluation when update failed
      }  # D1 > 0 test
    } # end if accpoint
    else { # no acceptable point
      if (trace > 0) cat("No acceptable point\n")
      if ( (ig == ilast) || (abs(gradproj) < (1 + abs(fmin))*ctrl$eps*ctrl$eps ) ) { # remove ig > 2
        # we reset to gradient and did new linesearch
        keepgoing <- FALSE  # no progress possible
        if (conv < 0) { # conv == -1 is used to indicate it is not set
          conv <- 0
        }
        msg <- "Converged"
        if (trace > 0) cat(msg, "\n")
      } # end ig == ilast
      else {
        ilast <- ig  # reset to gradient search
        if (trace > 0) cat("Reset to gradient search\n")
      }  # end else ig != ilast
    }  # end else no accpoint
  }  # end main loop  (while keepgoing)
  if (maximize) fmin <- (-1) * fmin
  if (trace > 0) cat("Seem to be done Rvmminb\n")
  msg <- "Rvmminb appears to have converged"
  counts <- c(ifn, ig)
  names(counts) <- c("function", "gradient")
  ans <- list(par, fmin, counts, convergence=conv, msg, bdmsk)
  names(ans) <- c("par", "value", "counts", "convergence", 
                  "message", "bdmsk")
  ans    #return(ans)
}  ## end of Rvmminb

grchk <- function(xpar, ffn, ggr, trace=0, testtol=(.Machine$double.eps)^(1/3), ...) {
  # check gradient code in ggr for function in ffn
  # assume function ffn already passes checks in fnchk
  if (is.null(trace)) trace <- 0 else trace<-trace
  if (is.null(testtol)) testtol<-(.Machine$double.eps)^(1/3) else testtol<- testtol
  gname <- deparse(substitute(ggr)) # ?? May be embedded in ugr
  if (trace>1) cat("Gradient test with tolerance = ",testtol,"\n")
  if (trace>1) cat("Analytic gradient uses function ",gname,"\n")
  fval <- ffn(xpar, ...)
  if (trace>1) {
    cat("function at parameters = ", fval," with attributes:\n")
    print(attributes(fval))
  }
  if (trace>1) cat("Compute analytic gradient\n")
  ga<-ggr(xpar, ...)
  if (trace>1) print(ga)
  if (trace>1) cat("Compute numeric gradient\n")
  gn <- gradient(fn=ffn, x=xpar, ...) # numerically approximated gradient
  if (trace>1) print(gn)
  # Now test for equality (090612: ?? There may be better choices for the tolerances.
  if (trace>0) {
    cat("gradient test tolerance = ",testtol,"  fval=",fval,"\n")
    cat(" compare to max(abs(gn-ga))/(1+abs(fval)) = ",max(abs(gn-ga))/(1+abs(fval)),"\n")
  }
  if (max(abs(gn-ga))/(1 + abs(fval)) >= testtol) { 
    if(trace>0) cat("Gradient function might be wrong - check it! \n")
    gradchk<-FALSE
  } else gradchk<-TRUE
  attr(gradchk, "ga")<-ga
  attr(gradchk, "gn")<-gn
  gradchk
} # end gradchk

bmchk <- function(par, lower = NULL, upper = NULL, 
                  bdmsk = NULL, trace = 0, tol = NULL, shift2bound = TRUE) {
  ## Bounds and masks check
  #  ?? check use of par and bvec -- can we simplify?
  # 20101031 -- issue of bounds not working correctly
  #  - -Inf seems to upset bounds
  #  - no upper bounds gives troubles (applies to Rcgmin too!)
  #
  # Input:
  #  par = a vector containing the starting point
  #  lower = vector of lower bounds on parameters
  #  upper = vector of upper bounds on parameters
  # Note: free parameters outside bounds will be adjusted to
  #   bounds.
  # bdmsk = control vector for bounds and masks. Parameters
  #   for which bdmsk are 1
  # are unconstrained or 'free', those with bdmsk 0 are
  #   masked i.e., fixed.
  # For historical reasons, we use the same array as an
  #   indicator that a parameter is at a lower bound (-3)
  #   or upper bound (-1)
  # trace = control of output: 0 for none (default), >0 for
  #   output
  # tol = tolerance for closeness of upper and lower bounds
  # shift2bound = TRUE if we adjust par values so they are
  #   feasible
  ##
  # Output:
  #    A list with components:
  #     bvec: The parameters adjusted to the nearest bound.
  #     bdmsk: adjusted input masks
  #     bchar: indicator for humans -- "-","L","F","U","+","M"
  #        for out-of-bounds-low, lower bound, free, 
  #            upper bound, out-of-bounds-high, masked (fixed)
  #     lower: adjusted lower bounds
  #     upper: adjusted upper bounds
  #     nolower: TRUE if no lower bounds, FALSE otherwise
  #     noupper: TRUE if no upper bounds, FALSE otherwise
  #     bounds:  TRUE if any bounds, FALSE otherwise
  #     admissible: TRUE if admissible, FALSE if not
  #        No lower bound exceeds an upper bound. That is the 
  #        bounds themselves are sensible. This condition has 
  #        nothing to do with the starting parameters.
  #     maskadded: TRUE if a mask is added, FALSE if not
  #        This implies very close upper and lower bounds for 
  #        the parameters. See the code for the implementation.
  #     parchanged: TRUE if parameters changed, FALSE if not
  #        parchanged = TRUE means that parameters are 
  #        INFEASIBLE, or they would not be changed.
  #     onbound: TRUE if any parameter equal to a bound
  #     feasible: TRUE if feasible, FALSE otherwise
  #
  ########## length of vectors #########
  n <- length(par)
  bvec <- par
  ############# bounds and masks ################
  # set default masks if not defined
  bchar <- rep("F",n) # make sure these are defined
  if (is.null(bdmsk)) {
    #        cat("setting bdmsk\n")
    bdmsk <- rep(1, n)
  }
  if (trace > 2) {
    cat("bdmsk:")
    print(bdmsk)
  }
  # check if there are bounds
  if (is.null(lower) || !any(is.finite(lower))) 
    nolower <- TRUE
  else nolower <- FALSE
  if (is.null(upper) || !any(is.finite(upper))) 
    noupper <- TRUE
  else noupper <- FALSE
  if (nolower && noupper && all(bdmsk == 1)) 
    bounds <- FALSE
  else bounds <- TRUE
  if (trace > 2) 
    cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, 
        " bounds = ", bounds, "\n")
  if (nolower) 
    lower <- rep(-Inf, n)
  if (noupper) 
    upper <- rep(Inf, n)
  
  ## adjust tolerance for masks and parameters ON bounds
  if (is.null(tol) || tol <= 0.0) {
    tol <- .Machine$double.eps * max(abs(par), 1) # use par as bounds Inf
  }
  ######## check bounds and masks #############
  parchanged <- FALSE  ## must be set BEFORE if (bounds) ...; 
  ## equivalent to feasible<-TRUE
  feasible <- TRUE
  admissible <- TRUE  # similarly must set before we look at bounds
  maskadded <- FALSE  # similarly set here
  onbound <- FALSE
  if (bounds) {
    # Make sure to expand lower and upper
    if (!nolower & (length(lower) < n)) 
    {
      if (length(lower) == 1) {
        lower <- rep(lower, n)
      }
      else {
        stop("1<length(lower)<n")
      }
    }  # else lower OK
    if (!noupper & (length(upper) < n)) 
    {
      if (length(upper) == 1) {
        upper <- rep(upper, n)
      }
      else {
        stop("1<length(upper)<n")
      }
    }  # else upper OK
    # At this point, we have full bounds in play
    ######## check admissibility ########
    if (any(lower[which(bdmsk != 0)] > upper[which(bdmsk != 0)])) admissible <- FALSE
    if (trace > 0) cat("admissible = ", admissible, "\n")
    if (any((upper - lower) < tol)) { # essentially masked
      makemask<-which(upper - lower < tol)
      warning("Masks (fixed parameters) set by bmchk due to tight bounds. CAUTION!!")
      if (trace > 0) {
        cat("Imposing mask as lower ~= upper for following parameters\n")
        print(makemask)
      }
      bdmsk[makemask] <- 0
      bchar[makemask] <- "M"
      # lower[makemask]<- -Inf
      # upper[makemask]<-  Inf
      maskadded <- TRUE
    }
    if (trace > 0) cat("maskadded = ", maskadded, "\n")
    ######## check feasibility ########
    if (admissible) {# This implementation as a loop, but try later to vectorize
      for (i in 1:n) {
        if (bdmsk[i] == 0) {
          bchar[i] <- "M"
          # NOTE: we do not change masked parameters, even if out of
          #   bounds
          if (!nolower) {
            if (bvec[i] < lower[i]) {
              if (trace > 0) {
                cat("WARNING: ", bvec[i], " = MASKED x[", 
                    i, "] < lower bound = ", lower[i], "\n")
              }
              feasible <- FALSE
            }
          }
          if (!noupper) {
            if (bvec[i] > upper[i]) {
              if (trace > 0){
                cat("WARNING: ", bvec[i], " = MASKED x[", 
                    i, "] > upper bound = ", upper[i], "\n")
              }
              feasible<-FALSE
            }
          }
        }
        else { # not masked, so must be free or active constraint
          if (!nolower) {
            if (bvec[i] <= lower[i]) {
              # Gave trouble 130924 -- <= not < -- bmtest in nlpor
              # changed 090814 to ensure bdmsk is set; 110105 < not <=
              if (bvec[i] != lower[i]){
                bdmsk[i] <- -3.5 # OUT OF BOUNDS LOW
                bchar[i] <- "-"
                feasible<-FALSE
                if (shift2bound) {
                  parchanged <- TRUE
                  if (trace > 0) cat("WARNING: x[", i, "], set ", 
                                     bvec[i], " to lower bound = ", lower[i], "\n")
                  bvec[i] <- lower[i]
                  bdmsk[i] <- -3 
                  bchar[i] <- "L"
                }
              } else {
                bdmsk[i] <- -3  # active lower bound
                bchar[i] <- "L"
              }
            }
          } # nolower
          if (!noupper) {
            if (bvec[i] >= upper[i]) {
              # changed 090814 to ensure bdmsk is set; 110105 > not >=
              if (bvec[i] != upper[i]){
                bdmsk[i] <- -0.5 # OUT OF BOUNDS HIGH
                bchar[i] <- "+"
                feasible<-FALSE
                if (shift2bound) {
                  parchanged <- TRUE
                  if (trace > 0) cat("WARNING: x[", i, "], set ",
                                     bvec[i], " to upper bound = ", upper[i], "\n")
                  bvec[i] <- upper[i]
                  bdmsk[i] <- -1   # active upper bound
                  bchar[i] <- "U"
                }
              } else { 
                bdmsk[i] <- -1  # active upper bound
                bchar[i] <- "U"
              }
            }
          } # noupper
        }  # end not masked
        ## Should NOT proceed with optimization if feasible and parchanged 
        ## both FALSE (out of bounds)
      }  # end loop for bound/mask check
    } # if admissible
  } # if bounds
  if (trace > 0) 
    cat("parchanged = ", parchanged, "\n")
  if (any(bvec == lower) || any(bvec == upper)){
    onbound <- TRUE
    if (trace > 0) cat("At least one parameter is on a bound\n")
  }
  ############## end bounds check #############
  bcout <- list(bvec, bdmsk, bchar, lower, upper, nolower, noupper, 
                bounds, admissible, maskadded, parchanged, feasible, onbound)
  names(bcout) <- c("bvec", "bdmsk", "bchar", "lower", "upper", "nolower", "noupper", 
                    "bounds", "admissible", "maskadded", "parchanged", "feasible", "onbound")
  # Note bdmsk, lower and upper are returned because they are
  #   modified (length, etc.)
  return(bcout)
}  ## end of bmchk.R 

