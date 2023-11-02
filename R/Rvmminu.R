Rvmminu <- function(par, fn, gr=NULL, control = list(), ...) {
  #
  #
  #  Author:  John C Nash
  #  Date:    Dec 6, 2021 update
  #
  ## An R version of the Nash version of Fletcher's Variable
  #   Metric minimization -- unconstrained parameters
  # This uses a simple backtracking line search.
  #
  # Input:
  #  par = a vector containing the starting point
  #  fn = objective function (assumed to be sufficeintly
  #    differentiable)
  # gr = gradient of objective function, provided as a function
  #   or the character name of a numerical approximation function
  #  
  # 
  #    This space to match Rvmminb comment line spacing
  # 
  # 
  # 
  #
  #
  #
  #
  # control = list of control parameters
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
  # message: A character string giving any additional
  #   information returned by the optimizer, or 'NULL'.
  #
  #
  #
  #
  #
  #
  #
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
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##################################################################
## Set working parameters (See CNM Alg 22)
  if (trace > 0) 
     cat("Rvmminu -- J C Nash 2009-2015 - an R implementation of Alg 21\n")
  bvec <- par  # copy the parameter vector
  n <- length(bvec)  # number of elements in par vector
  if (trace > 0) {
     cat("Problem of size n=", n, "  Dot arguments:\n")
     print(fargs)
  }
  ifn <- 1  # count function evaluations
#  stepredn <- 0.2  # Step reduction in line search
#  reltest <- 100  # relative equality test
  ceps <- .Machine$double.eps * reltest
  dblmax <- .Machine$double.xmax  # used to flag bad function
  #############################################
  # gr MUST be provided
  if (is.null(gr)) {  # if gr function is not provided STOP 
    stop("A gradient calculation (analytic or numerical) MUST be provided for Rvmminu") 
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
  ############# end setup gr ####################
  #
  f<-try(fn(bvec, ...), silent=FALSE) # Compute the function.
  if (inherits(f, "try-error") | is.na(f) | is.null(f) | is.infinite(f)) {
     msg <- "Initial point gives inadmissible function value"
     conv <- 20
     if (trace > 0) 
        cat(msg, "\n") # change NA to dblmax 110524
     ans <- list(bvec, dblmax, c(ifn, 0), conv, msg)  #
     names(ans) <- c("par", "value", "counts", "convergence", 
       "message")
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
#
#
#
#  vertical spacing for bounds and masks adjustment of gradient 
#  in Rvmminb
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
      t <- as.vector(-B %*% g)  # compute search direction
      if (!all(is.numeric(t))) 
          t <- rep(0, n)  # 110619
      if (trace > 2) {
          cat("t:")
          print(t)
      }
#  space for masks and box constraints        
#      if (trace > 2) {
#          cat("adj-t:")
#          print(t)
#      }
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
        while (changed && (!accpoint)) {
          # We seek a lower point, but must change parameters too
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
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
            accpoint <- (f < fmin + gradproj * steplength * acctol) # NOTE: < not <=
            if (trace > 2) cat("accpoint = ", accpoint,"\n")
            if (! accpoint) {
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
#
#
#
#  Rvmmminb uses this vertical space for reactivating constraints
#
#
#
#
#
#
#
#
#
#
#
#
#
        test <- try(g <- mygr(bvec, ...), silent = FALSE)
        if (inherits(test, "try-error") ) stop("Bad gradient!!")
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
# for Rvmminb active mask or constraint adjust g
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
    if (trace > 0) cat("Seem to be done Rvmminu\n")
    msg <- "Rvmminu appears to have converged"
    counts <- c(ifn, ig)
    names(counts) <- c("function", "gradient")
    ans <- list(par, fmin, counts, convergence=conv, msg)
    names(ans) <- c("par", "value", "counts", "convergence", 
        "message")
    ans    #return(ans)
}  ## end of Rvmminu
