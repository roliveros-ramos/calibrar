.Rvmminbxx = function(par, fn, gr = NULL, lower = NULL, 
                      upper = NULL, bdmsk = NULL, control = list(), ...) {
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
  if ((inherits(f, "try-error")) | is.na(f) | is.null(f) | is.infinite(f)) {
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
  
  BDMSK = matrix(nrow=maxit, ncol=n)
  jj = 1
  
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
  
  fmin <- f # ROR: new init
  
  while (keepgoing) { ## main loop -- must remember to break out of it!
    if (ilast == ig) { # reset the approx. inverse hessian B to unit matrix
      B <- diag(1, n, n)  # create unit matrix of order n
      if (trace > 1) cat("Reset Inv. Hessian approx at ilast = ", ilast, "\n")
    }
    # fmin <- f
    if (trace > 0) cat(" ", ifn, " ", ig, " ", fmin, "\n")
    # par <- bvec  # save parameters
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
          if (inherits(f, "try-error")) f <- .Machine$double.xmax
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
      if (inherits(test, "try-error")) stop("Bad gradient!!")
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
    
    BDMSK[jj, ] = bdmsk
    jj = jj + 1
    
  }  # end main loop  (while keepgoing)
  
  BDMSK = BDMSK[1:(jj-1), ]
  
  if (maximize) fmin <- (-1) * fmin
  if (trace > 0) cat("Seem to be done Rvmminb\n")
  msg <- "Rvmminb appears to have converged"
  counts <- c(ifn, ig)
  names(counts) <- c("function", "gradient")
  ans <- list(par, fmin, counts, convergence=conv, msg, bdmsk, trace=list(bdmsk=BDMSK))
  names(ans) <- c("par", "value", "counts", "convergence", 
                  "message", "bdmsk", "trace")
  ans    #return(ans)
}  ## end of Rvmminb
