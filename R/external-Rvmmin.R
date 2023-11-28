Rvmmin <- function(par, fn, gr = NULL, lower = NULL, 
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
  #    checkgrad = FALSE by default. Check analytic gradient 
  #            against numDeriv results.
  #    checkbounds = TRUE by default. Check parameters and bounds
  #            for addmissible bounds and feasible start.
  #    keepinputpar = FALSE if Rvmmin will move out-of-bound parameters
  #            to nearest bound WITH WARNING. This is default behaviour.
  #            If TRUE, stop execution. That is, do not allow starts to
  #            be modified by bmchk function. This is only in Rvmmin(), not
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
        if (is.null(control$checkgrad)) control$checkgrad <- FALSE
        if (control$checkgrad) { # check gradient
           testgrad<-grchk(par, fn, gr, trace=control$trace, ...)
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
    ans <- Rvmminb(par, fn, gr, lower = lower, 
        upper = upper, bdmsk = bdmsk, control = control, ...)
    } else {
       ans <- Rvmminu(par, fn, gr, control = control, ...)
    } #   return(ans) 
}  ## end of Rvmmin
