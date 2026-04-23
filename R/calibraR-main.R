# calibrar package: Automated Calibration for Complex Models --------

#' calibrar package: Automated Calibration for Complex Models
#'
#'\if{html}{\figure{logo_small.png}{options: style='float: right' alt='logo' width='150'}}
#'
#' This package allows the parameter estimation (i.e. calibration) of complex models, including stochastic ones. It implements generic functions that can be used for fitting any type of models, especially those with non-differentiable objective functions, with the same syntax as base::optim. 
#' It supports multiple phases estimation (sequential parameter masking), constrained optimization (bounding box restrictions) and automatic parallel computation of numerical gradients. 
#' Some common maximum likelihood estimation methods and automated construction of the objective function from simulated model outputs is provided.  
#' See <https://roliveros-ramos.github.io/calibrar/> for more details.
#' 
#' @title Automated Calibration for Complex Models
#' @description Automated Calibration for Complex Models
#' @name calibrar-package
#' @aliases calibrar-package calibrar
#' @author Ricardo Oliveros-Ramos Maintainer: Ricardo Oliveros-Ramos
#' <ricardo.oliveros@@gmail.com>
#' @references calibrar: an R package for the calibration of ecological models (Oliveros-Ramos and Shin 2014)
#' @keywords calibration
#' @examples
#' \dontrun{
#' require(calibrar)
#' set.seed(880820)
#' path = NULL # NULL to use the current directory
#' # create the demonstration files
#' demo = calibrar_demo(model="PoissonMixedModel", L=5, T=100) 
#' # get calibration information
#' calibrationInfo = calibration_setup(file=demo$path)
#' # get observed data
#' observed = calibration_data(setup=calibrationInfo, path=demo$path)
#' # read forcings for the model
#' forcing = read.csv(file.path(demo$path, "master", "environment.csv"), row.names=1)
#' # Defining 'runModel' function
#' runModel = function(par, forcing) {
#' output = calibrar:::.PoissonMixedModel(par=par, forcing=forcing)
#' # adding gamma parameters for penalties
#' output = c(output, list(gammas=par$gamma)) 
#' return(output)
#' }
#' # real parameters
#' cat("Real parameters used to simulate data\n")
#' print(demo$par)
#' # objective functions
#' obj  = calibration_objFn(model=runModel, setup=calibrationInfo, 
#'                                observed=observed, forcing=forcing)
#' cat("Starting calibration...\n")
#' control = list(weights=calibrationInfo$weights, maxit=3.6e5) # control parameters
#' cat("Running optimization algorithms\n", "\t", date(), "\n")
#' cat("Running optim AHR-ES\n")
#' ahr = calibrate(par=demo$guess, fn=obj, lower=demo$lower, upper=demo$upper, control=control)
#' summary(ahr)
#' } 
"_PACKAGE"

# calibrate ---------------------------------------------------------------

#' Sequential parameter estimation for the calibration of complex models
#'
#' \code{calibrate()} minimises an objective function \code{fn} for the calibration of
#' complex (possibly expensive and stochastic) models. It supports sequential calibration
#' in multiple phases (progressively activating parameters), replicated evaluations for
#' stochastic objectives, restartable runs, and parallel execution patterns for
#' computationally intensive models.
#'
#' @details
#' \strong{Sequential phases.} When \code{phases} is provided, parameters are activated
#' progressively across phases. In phase \eqn{k}, parameters with \code{phases <= k} are
#' active and estimated, while all others are held fixed. If \code{phases} is omitted, all
#' parameters are estimated in a single phase.
#'
#' \strong{Replicated evaluations.} The argument \code{replicates} controls the number of
#' replicate evaluations of \code{fn} per phase. It can be a single integer (applied to all
#' phases) or a vector of length equal to the number of phases. Replication is useful for
#' stochastic models to reduce Monte Carlo noise and/or to target robust parameter sets.
#'
#' \strong{Default method selection.} If \code{method} is not provided, \code{calibrate()}
#' defaults to \code{"Rvmmin"} when \code{all(replicates == 1)} (deterministic objective),
#' and to \code{"AHR-ES"} otherwise (replicated/stochastic objective).
#'
#' \strong{Multi-objective outputs.} The objective function \code{fn} may return a scalar
#' (single-objective) or a numeric vector of objective components (multi-objective).
#' Multi-objective optimisation is currently supported only by methods in
#' \code{multiMethods} (e.g., \code{"AHR-ES"}). For single-objective methods, \code{fn} must be
#' scalar; alternatively, objective components can be aggregated into a scalar by defining
#' an aggregated objective in \code{calibration_objFn(aggregate = TRUE)}.
#'
#' \strong{Aggregation.} When \code{fn} is created via \code{calibration_objFn()}, metadata
#' such as the number of components (\code{nvar}) and component weights (\code{weights}) can
#' be stored as attributes. If \code{aggregate = TRUE}, \code{fn} is scalar and can be used
#' with any method. If \code{aggregate = FALSE} (vector output), only multi-objective methods
#' can use it directly.
#'
#' \strong{Parallel execution and run directories.} If \code{parallel = TRUE},
#' \code{calibrate()} may distribute replicate evaluations and/or finite-difference
#' gradient computations across multiple cores. The number of cores is controlled by
#' \code{control$ncores}. For file-based or externally executed models, \code{control$master}
#' and \code{control$run} can be used to manage a master/template directory and per-run
#' working directories (created if needed).
#'
#' \strong{Restart.} Partial results can be written and used to resume a calibration via
#' \code{control$restart.file}. Restart functionality is currently available only for
#' \code{"AHR-ES"}, \code{"Rvmmin"}, and \code{"hjn"}.
#'
#' @section Choosing an optimisation method:
#' For smooth deterministic objectives, gradient-based methods such as \code{"Rvmmin"},
#' \code{"L-BFGS-B"}, \code{"LBFGSB3"}, or \code{"nlminb"} are often efficient (with box
#' constraints when needed). For noisy or stochastic objectives (typically when
#' \code{replicates > 1}), heuristic/global methods such as \code{"AHR-ES"} are generally
#' more appropriate. Not all methods are directly comparable (local vs.\ global,
#' deterministic vs.\ stochastic, scalar vs.\ vector-valued objectives); method choice
#' should reflect the structure of \code{fn} and the presence of constraints.
#'
#' @section Notes:
#' \code{"SANN"} is included for compatibility with \code{stats::optim()}, but it is highly
#' sensitive to tuning and often performs poorly on continuous problems under default
#' settings. For stochastic objectives, \code{"AHR-ES"} is the intended default within
#' \code{calibrate()}.
#'
#' @param par A numeric vector or list. The length of \code{par} defines the number of
#' parameters to be estimated (i.e., the dimension of the problem).
#' @param fn The objective function to be minimised. It should accept a parameter vector
#' (or list, depending on the wrapper) as first argument and return either a scalar value
#' (single-objective) or a numeric vector (multi-objective).
#' @param gr A function computing the gradient of \code{fn}. If \code{NULL}, a numerical
#' approximation is used. Alternatively, a character string can specify the numerical
#' gradient scheme: \code{"central"}, \code{"forward"} (default), \code{"backward"}, or
#' \code{"richardson"}.
#' @param ... Additional arguments passed to \code{fn} and \code{gr}.
#' @param method Optimisation method(s) to be used. Can be a single method name or a
#' vector of method names (e.g., one per phase). If \code{NULL}, a default is chosen based
#' on \code{replicates} (see Details).
#' @param lower Lower bounds for parameters. One value or a vector of the same length as
#' \code{par}. \code{NA} is treated as \code{-Inf}. Default is unconstrained.
#' @param upper Upper bounds for parameters. One value or a vector of the same length as
#' \code{par}. \code{NA} is treated as \code{Inf}. Default is unconstrained.
#' @param phases Optional integer vector of the same length as \code{par}, indicating the
#' phase at which each parameter becomes active. If omitted, all parameters are active in
#' a single phase.
#' @param control A list of control options. Common options include \code{ncores},
#' \code{run}, \code{master}, \code{verbose}, \code{REPORT}, \code{restart.file},
#' \code{gradient}, and \code{gr.method}. Additional solver-specific options may be passed
#' through to the underlying optimiser.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Currently not implemented.
#' @param replicates Integer or integer vector controlling the number of replicate
#' evaluations of \code{fn} per phase. The default is \code{1} (deterministic objective).
#' @param parallel Logical. Enable parallel computation (e.g., for replicated evaluations
#' and/or numerical gradients) using up to \code{control$ncores} cores.
#'
#' @return
#' An object of class \code{"calibrar.results"} with components:
#' \describe{
#' \item{par}{Best parameter values found (returned with the same structure as input \code{par}).}
#' \item{value}{Objective value at \code{par}.}
#' \item{counts}{Number of calls to \code{fn} and \code{gr} (where applicable).}
#' \item{convergence}{Convergence code returned by the underlying optimiser.}
#' \item{message}{Additional information returned by the optimiser, if any.}
#' \item{method}{The optimisation method used in the final phase.}
#' \item{fn}{The objective function.}
#' \item{active}{Logical vector indicating which parameters were active (estimated).}
#' \item{elapsed}{Elapsed time for the full calibration run.}
#' \item{trace}{Tracing information and per-phase outputs (when available).}
#' }
#'
#' @seealso \code{\link{calibration_setup}}, \code{\link{calibration_data}},
#' \code{\link{calibration_objFn}}
#'
#' @examples
#' calibrate(par=rep(NA, 5), fn=sphereN)
#' \dontrun{
#' calibrate(par=rep(NA, 5), fn=sphereN, replicates=3)
#' calibrate(par=rep(0.5, 5), fn=sphereN, replicates=3, lower=-5, upper=5)
#' calibrate(par=rep(0.5, 5), fn=sphereN, replicates=3, lower=-5, upper=5, phases=c(1,1,1,2,3))
#' calibrate(par=rep(0.5, 5), fn=sphereN, replicates=c(1,1,4), lower=-5, upper=5, phases=c(1,1,1,2,3))
#' }
#'
#' @author Ricardo Oliveros-Ramos
#' @family optimisers
#' @export
calibrate = function(par, fn, gr, ..., method, lower, upper, phases, control, 
                     hessian, replicates, parallel) {
  obj = if (!missing(fn)) fn else par
  UseMethod("calibrate", obj)
}

#' @export
#' @rdname calibrate
calibrate.default = function(par, fn, gr = NULL, ..., method = NULL,  
                             lower = NULL, upper = NULL, phases = NULL, 
                             control = list(), 
                             hessian = FALSE, replicates=1, parallel=FALSE) {
  
  if(is.null(method)) {
    method = if(all(replicates==1)) "Rvmmin" else "AHR-ES"
  }
  
  methods = c("L-BFGS-B", "nlminb", "Rcgmin", "Rvmmin", 
              "hjn", "spg", "LBFGSB3", "CMA-ES", "genSA", "DE", "soma", "genoud", "PSO", 
              "hybridPSO", "mads", "hjkb", "nmkb", "bobyqa", "AHR-ES", "Nelder-Mead", 
              "CG", "BFGS", "SANN")
  
  method = match.arg(method, choices=methods, several.ok = TRUE)
  
  # calibrate should dispatch on 'fn'
  if(inherits(fn, "objFn")) {
    control$nvar = attr(fn, "nvar")
    control$weights = attr(fn, "weights")
    names(control$weights) = attr(fn, "variables")
  }
  
  # check function and method
  
  if(!all(method %in% multiMethods)) {
  
    if(inherits(fn, "objFn")) {
      agg = attr(fn, "aggregate")
      if(is.null(agg)) warning("Update your objective function to the last version of the 'calibrar' package.")
      if(!isTRUE(agg)) 
        stop(sprintf("Method '%s' does not support multi-objective optimization, use aggregate=TRUE in 'calibration_objFn'.",
                     method))
    }
    
    if(!is.null(control$nvar)) 
      if(control$nvar != 1) stop(sprintf("Method '%s' does not support multi-objective optimization, 'fn' must be scalar.",
                                         method))
    control$nvar = 1
    
  }
  
  fnx = NULL
  if(inherits(fn, "objFn")) {
    control$aggFn = attr(fn, "aggFn")
    fnx = attr(fn, "fn")
    if(!is.null(fnx)) fnx = match.fun(fnx)
  }

  # if 'gr' is a character, assume it specify the type of numerical approximation
  if(!is.null(gr) & is.character(gr)) {
    if(!is.null(control$gr.method))  
      warning(sprintf("Ignoring provided gr.method='%s', using gr='%s'.",
                      control$gr.method, gr))
    control$gr.method = gr
    gr = NULL
  }
  
  control = .checkControl_calibrate(control=control, method=method, par=par, fn=fn, ...)
  # here, we calculate 'batchsize' according to the method, copyMaster is called just once.
  
  # check for a partial results file to restart from a completed phase
  restart = .restartCalibration(control, type="partial")
  
  skeleton = as.relistable(par)
  par    = unlist(par)
  lower  = unlist(lower)
  upper  = unlist(upper)
  phases = unlist(phases)
  
  npar = length(par)
  
  # checking conformance of all arguments
  if(is.null(names(par))) names(par) = .printSeq(npar, preffix="par")
  
  fn = match.fun(fn)
  if(!is.null(gr)) gr = match.fun(gr)

  phases     = .checkPhases(phases=phases, npar=npar)
  bounds     = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess      = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  
  par     = guess
  lower   = bounds$lower
  upper   = bounds$upper
  nphases = max(phases, na.rm=TRUE)
  
  replicates = .checkReplicates(replicates, nphases)
  control$repmax = max(replicates, na.rm=TRUE)
  
  output = if(isTRUE(restart)) .getResults(control=control, type="partial") else list(phase=0)
  output$phase = output$phase + 1
 
  msg0 = sprintf("\nStarting calibration from phase %d.\n", output$phase)
  msg1 = sprintf("\nRestarting calibration from phase %d.\n", output$phase)

  if(isTRUE(restart)) message(msg1)
  if(!isTRUE(restart) & isTRUE(control$verbose)) message(msg0)
  
  conv = .checkConvergence(control, nphases, method)
 
  control$parallel = parallel
  
  if(!is.null(control$master)) message(sprintf("Using 'master' directory: %s", control$master))
  if(!is.null(control$run)) message(sprintf("Using 'run' directory: %s", control$run))
   
  xtm1 = Sys.time() # timer for all phases
  
  # -------------------------------------------
  # start the sequential parameter estimation
  # -------------------------------------------
  if(nphases > 1) 
    # message(sprintf("PARAMETER ESTIMATION IN %s PHASES.", toupper(make_number(nphases))))
    message(sprintf("Parameter estimation in %s phases.", make_number(nphases)))
  
  for(phase in seq(from=output$phase, to=nphases)) {
    
    if(output$phase > nphases) break
    
    control$maxgen      = conv$maxgen[phase]
    control$maxit       = conv$maxit[phase]
    control$convergence = conv$convergence[phase]
    control$abstol      = conv$abstol[phase]
    control$reltol      = conv$reltol[phase]
    method              = conv$method[phase]
    
    active = (phases <= phase) # NAs are corrected in .calibrar 
    
    msg2 = sprintf("\n- Phase %d: %d out of %d parameters are currently active.", 
                   phase, sum(active, na.rm=TRUE), npar)
    msg3 = sprintf("\tStarted at %s.\n", date())
    if(nphases > 1) message(msg2)
    if(isTRUE(control$verbose)) message(msg3)
    if(nphases==1) message(sprintf("Using optimization method '%s'.", method))
    if(nphases!=1) message(sprintf("\tUsing optimization method '%s'.", method))
    
    # call optimisers using .calibrar handler
    tm1 = Sys.time()
    temp = .calibrar(par=par, fn=fn, gr = gr, ..., method = method, 
                   lower = lower, upper = upper, control = control, 
                   hessian = hessian, active=active, skeleton=skeleton,
                   replicates=replicates[phase])
    tm2 = Sys.time()
    
    output$phases[[phase]] = temp # trim?
    output$phase = phase
    
    .createOutputFile(output, control, type="partial", phase=phase) 
    
    par = temp$par #update parameter guess
    control = .updateControl(control=control, opt=temp, method=method)  # update CVs? 
    
    if(!is.null(control$restart.file) & !is.null(fnx)) {
      
      wd = getwd()
      xpar = relist(flesh = par, skeleton = skeleton)
      ifile = sprintf(".%s-phase%d.par", control$restart.file, phase)
      saveRDS(xpar, file=ifile)
      .setWorkDir(control$run, i=0)
      partial = try(fnx(xpar))
      setwd(wd)
      if(!inherits(partial, "try-error")) {
        ifile = sprintf("%s-phase%d.simulated", control$restart.file, phase)
        saveRDS(partial, file=ifile)
      }
      
    } 
    
    xmsg = temp$message
    if(is.null(xmsg)) xmsg = NA
    if(is.na(xmsg)) xmsg = "-"
    
    msg = c(sprintf("Phase %d finished (%s)", phase, format_difftime(tm1, tm2)),
            sprintf("Elapsed time: %s", format_difftime(tm1, tm2)),
            sprintf("Function value: %g", temp$value),
            paste(c("Parameter values:",sprintf("%0.3g", par[which(active)])), collapse=" "))
    
    msg = if(nphases == 1) msg[-1] else paste("\t", msg[-2], sep="")
    msg = if(phase == nphases) c(msg, sprintf("\nStatus: %s", xmsg)) else msg
    msg = paste(msg, collapse="\n")
    message(msg)
    
  } # end of phases
  xtm2 = Sys.time() # timer for all phases
  
  isActive = (phases>0) & !is.na(phases)
  paropt = output$phases[[nphases]]$par # save parameters of last phase
  
  paropt = relist(paropt, skeleton)
  class(paropt) = setdiff(class(paropt), "relistable")
  
  trace = output$phases[[nphases]]$trace # must be a list
  if(is.null(trace)) trace = list()
  if(!is.null(output$phases[[nphases]]$generations))
    trace$generations = output$phases[[nphases]]$generations
  trace = c(trace, list(partial=output$phases[[nphases]]$partial), output)
  
  final = list(par         = paropt, 
               value       = output$phases[[nphases]]$value, 
               counts      = output$phases[[nphases]]$counts,
               convergence = output$phases[[nphases]]$convergence,
               message     = output$phases[[nphases]]$message,
               method      = method,
               fn          = fn, 
               active      = isActive,
               elapsed     = format_difftime(xtm1, xtm2, value = TRUE),
               trace       = trace)
  
  # output = c(final, output)
  class(final) = c("calibrar.results")
  .createOutputFile(final, control, type="results")
  
  return(final)
  
}


#' @export
calibrate.TMB = function(par, fn, gr = NULL, ..., method = NULL,  
                         lower = NULL, upper = NULL, phases = NULL, 
                         control = list(), 
                         hessian = NULL, replicates=1, parallel=FALSE) {
  
  methods = c("nlminb", "BFGS", "L-BFGS-B", "Rcgmin", "Rvmmin", "spg", 
              "LBFGSB3", "AHR-ES")
  
  if(!is.null(method)) fn$method = method
  
  method = match.arg(fn$method, choices=methods)
  
  if(!is.null(hessian)) fn$hessian = hessian
    
  out = calibrate(par=par, fn=fn$fn, gr=fn$gr, lower=lower, upper=upper, 
                  phases=phases, method=method, control=control, 
                  hessian=fn$hessian, replicates=1, parallel=parallel)
  
  return(out)
  
}

# optim2 ------------------------------------------------------------------


#' Unified optimisation interface with structured parameters and parallel numerical gradients
#'
#' \code{optim2()} provides a unified interface to multiple deterministic and stochastic
#' optimisation methods, combining the algorithms available through \code{stats::optim()}
#' with a small set of additional solvers accessed via \pkg{calibrar}'s internal dispatcher.
#'
#' @details
#' \strong{Methods.} The current selection includes
#' (i) base \code{stats::optim()} methods, (ii) additional gradient-based solvers from external
#' packages, and (iii) heuristic/global methods.
#'
#' \strong{Comparability.} Not all methods are directly comparable: some are deterministic
#' local optimisers (e.g., quasi-Newton), others are derivative-free local searches, and others
#' are stochastic/global heuristics. Method choice should reflect the objective function
#' (smooth vs. non-smooth, deterministic vs. stochastic) and the presence of constraints.
#'
#' \strong{Control arguments.} \code{optim2()} standardises a set of common control arguments
#' (e.g., iteration limits, tolerances, tracing, and scaling) where supported, while still
#' allowing method-specific control parameters to be passed through to the underlying solver.
#'
#' \strong{Parallel numerical gradients.} When analytic gradients are not provided,
#' \code{optim2()} can compute finite-difference numerical gradients and distribute function
#' evaluations across multiple cores, which can substantially reduce wall-clock time for
#' expensive objective functions.
#'
#' @section Choosing an optimisation method:
#' For smooth deterministic objectives, quasi-Newton methods (e.g., \code{"BFGS"} or
#' \code{"L-BFGS-B"} with box constraints) are often efficient when gradients are available or
#' can be reliably approximated. For box-constrained problems, consider \code{"L-BFGS-B"} or
#' the extended bounded solvers (e.g., \code{"Rvmmin"}, \code{"spg"}). For noisy or stochastic
#' objectives, heuristic/global methods (e.g., \code{"AHR-ES"}) may be more appropriate.
#'
#' @section Notes:
#' \code{"SANN"} is included for compatibility with \code{stats::optim()}, but it is highly
#' sensitive to tuning and often performs poorly on continuous problems under default
#' settings. For noisy or rugged objective functions, \code{"AHR-ES"} is generally the
#' recommended heuristic alternative within \code{optim2()}.
#'
#' @param active Boolean vector of the same length as \code{par}, indicating if the
#' parameter is used in the optimisation (\code{TRUE}) or held at a fixed value (\code{FALSE}).
#'
#' @return
#' A list with components:
#' \describe{
#' \item{par}{The best set of parameters found.}
#' \item{value}{The value of \code{fn} corresponding to \code{par}.}
#' \item{counts}{A two-element integer vector giving the number of calls to \code{fn} and
#' \code{gr} respectively. This excludes those calls needed to compute the Hessian, if
#' requested, and any calls to \code{fn} to compute a finite-difference approximation to the
#' gradient.}
#' \item{convergence}{An integer code. \code{0} indicates successful completion.}
#' \item{message}{A character string giving any additional information returned by the
#' optimizer, or \code{NULL}.}
#' \item{hessian}{Only if argument \code{hessian} is \code{TRUE}. A symmetric matrix giving an
#' estimate of the Hessian at the solution found. Note that this is the Hessian of the
#' unconstrained problem even if the box constraints are active.}
#' }
#'
#' @seealso \code{\link[stats]{optim}}, \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}}
#'
#' @examples
#' optim2(par=rep(NA, 5), fn=sphereN)
#' @inheritParams calibrate
#' @author Ricardo Oliveros-Ramos
#' @family optimisers
#' @export
optim2 = function(par, fn, gr = NULL, ..., 
                  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                             "Brent", "nlm", "nlminb", "Rcgmin", "Rvmmin", "hjn", 
                             "spg", "LBFGSB3", "AHR-ES"),
                  lower = -Inf, upper = +Inf, active = NULL, 
                  control = list(), hessian = FALSE, parallel=FALSE) {
  
  method = match.arg(method)
  
  out = .optim2(par=par, fn=fn, gr=gr, ..., lower=lower, upper=upper, active=active, 
                method=method, control=control, hessian=hessian, parallel=parallel)
  
  return(out)
    
}

# optimh ------------------------------------------------------------------

#' General-purpose optimisation with heuristic and derivative-free methods
#'
#' `optimh()` provides a unified interface to the heuristic and derivative-free
#' optimisation methods available through **calibrar**. It is intended for
#' objective functions that are expensive, noisy, non-smooth, stochastic, or
#' otherwise difficult to handle with standard gradient-based optimisers.
#'
#' @details
#' `optimh()` uses the same interface conventions as [optim2()] and
#' [calibrate()]: parameters may be supplied as a vector or list, box
#' constraints can be enforced through `lower` and `upper`, and subsets of
#' parameters can be held fixed through `active`.
#'
#' The available methods include both genuinely heuristic/global algorithms
#' (for example `"AHR-ES"`, `"CMA-ES"`, `"DE"`, `"genSA"`, `"PSO"`) and a few
#' derivative-free local searches (for example `"Nelder-Mead"`, `"bobyqa"`,
#' `"hjn"`, `"nmk"`). The function name should therefore be read broadly: it is
#' the package entry point for optimisation methods that are typically useful
#' when classical smooth, gradient-based optimisation is not the best fit.
#'
#' Analytic gradients are accepted for interface consistency, but many methods
#' implemented through `optimh()` ignore `gr`. Likewise, not all methods return
#' Hessian information even when `hessian = TRUE` is requested.
#'
#' @section Choosing a method:
#' Heuristic and evolutionary methods are often useful for rugged,
#' multi-modal, or stochastic objective functions. Derivative-free local
#' methods may be a better first choice when the objective is deterministic but
#' gradients are unavailable or unreliable. Method choice should still reflect
#' the structure of the problem: locality versus global exploration, box
#' constraints, stochasticity, and computational cost.
#'
#' @param gr Function computing the gradient of `fn`. It is included for API
#' consistency with [optim2()] and [calibrate()], but is ignored by most
#' methods available through `optimh()`.
#' @param method Optimisation method to use. Available choices are
#' `"AHR-ES"`, `"Nelder-Mead"`, `"SANN"`, `"hjn"`, `"bobyqa"`, `"CMA-ES"`,
#' `"genSA"`, `"DE"`, `"soma"`, `"genoud"`, `"PSO"`, `"hybridPSO"`, `"mads"`,
#' `"hjk"`, `"hjkb"`, `"nmk"`, and `"nmkb"`.
#'
#' @return
#' A list in the same general format as [optim2()], typically including:
#' \describe{
#' \item{par}{The best parameter values found.}
#' \item{value}{The value of `fn` at `par`.}
#' \item{counts}{Number of function and gradient evaluations reported by the
#' underlying solver, when available.}
#' \item{convergence}{Convergence code returned by the underlying solver.}
#' \item{message}{Additional information returned by the optimiser, if any.}
#' \item{hessian}{A Hessian estimate when supported and requested; otherwise
#' usually `NULL`.}
#' }
#'
#' @seealso [optim2()], [calibrate()], [ahres()], [gradient()]
#'
#' @examples
#' \dontrun{
#' optimh(
#'   par = rep(0.5, 5),
#'   fn = sphereN,
#'   lower = -5,
#'   upper = 5,
#'   method = "AHR-ES",
#'   control = list(maxit = 200)
#' )
#' }
#'
#' @inheritParams optim2
#' @author Ricardo Oliveros-Ramos
#' @family optimisers
#' @export
optimh = function(par, fn, gr = NULL, ...,  
                  method = c("AHR-ES", "Nelder-Mead", "SANN", "hjn", "bobyqa",
                             "CMA-ES", "genSA", "DE", "soma", "genoud", "PSO", 
                             "hybridPSO", "mads", "hjk", "hjkb", "nmk", "nmkb"), 
                  lower = -Inf, upper = +Inf, active = NULL,
                  control = list(), hessian = FALSE, parallel=FALSE) {
  
  method = match.arg(method)
  
  out = .optim2(par=par, fn=fn, gr=gr, ..., lower=lower, upper=upper, active=active, 
                method=method, control=control, hessian=hessian, parallel=parallel)
  
  return(out)
  
}

# ahres -------------------------------------------------------------------

#' Adaptive Hierarchical Recombination Evolutionary Strategy
#'
#' `ahres()` runs the AHR-ES optimiser directly through a dedicated interface.
#' It is primarily intended for bounded, derivative-free optimisation of
#' expensive or black-box objective functions, including calibration problems in
#' which the objective may be noisy or obtained from an external model.
#'
#' @details
#' `ahres()` accepts the same broad interface as the other optimisation entry
#' points in **calibrar**, but is specialised to the `"AHR-ES"` algorithm.
#' Parameters may be supplied as a vector or list. Internally, list-structured
#' parameters are flattened for optimisation and then restored before the result
#' is returned.
#'
#' The function supports box constraints through `lower` and `upper`, and
#' parameter masking through `active`, so only selected parameters are evolved
#' while the others are held fixed at their supplied starting values.
#'
#' If `control$parscale` and/or `control$fnscale` are supplied, optimisation is
#' performed on the corresponding scaled parameter space and objective scale,
#' and results are transformed back before being returned.
#'
#' `gr` is accepted only for interface compatibility and is ignored by this
#' method. In the current implementation, `control$maxit` or `control$maxgen`
#' must be supplied; otherwise the optimisation stops with an error.
#'
#' For scalar objectives, the default settings are usually sufficient. For
#' vector-valued objectives, advanced users may also need to provide suitable
#' `control` entries such as `weights`, `nvar`, and `aggFn`, depending on how
#' the objective has been constructed.
#'
#' @section Parallel execution:
#' AHR-ES parallelism is controlled through entries inside `control`, notably
#' `control$parallel` and `control$ncores`. The top-level `parallel` argument is
#' currently retained for interface consistency but is not used directly inside
#' `ahres()`.
#'
#' @param gr Function computing the gradient of `fn`. Ignored by `ahres()`.
#'
#' @return
#' An object of class `"ahres.results"`, typically with components:
#' \describe{
#' \item{par}{Best parameter values found, returned on the original scale.}
#' \item{value}{Objective value at `par`, returned on the original objective scale.}
#' \item{counts}{Number of function and gradient evaluations.}
#' \item{convergence}{Convergence code (`0` usually indicates a regular stop,
#' `1` that the maximum number of generations/evaluations was reached).}
#' \item{message}{A short message describing the stopping condition.}
#' \item{generations}{Number of generations performed.}
#' \item{partial}{Optional vector of objective components when the objective
#' function returns more than one value.}
#' \item{trace}{Optional tracing information when tracing is enabled.}
#' \item{active}{Information about active parameters used during optimisation.}
#' }
#'
#' @seealso [optimh()], [optim2()], [calibrate()], [summary.ahres.results()]
#'
#' @examples
#' \dontrun{
#' ahres(
#'   par = rep(0.5, 5),
#'   fn = sphereN,
#'   lower = -5,
#'   upper = 5,
#'   control = list(maxit = 200)
#' )
#' }
#'
#' @inheritParams optim2
#' @author Ricardo Oliveros-Ramos
#' @family optimisers
#' @export
ahres = function(par, fn, gr = NULL, ..., lower = -Inf, upper = +Inf, active = NULL, 
                 control = list(), hessian = FALSE, parallel=FALSE) {
  
  # par can be a list
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
  
  # Checking par
  active     = .checkActive(active=active, npar=npar)
  bounds     = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess      = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  
  par     = guess
  lower   = bounds$lower
  upper   = bounds$upper
  
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  if(is.null(names(par))) names(par) = .printSeq(npar, preffix="par")
  
  # getting functions
  fn = match.fun(fn)
  if(!is.null(gr)) warning("Gradient ignored by this method.")
  
  # when 'parscale' and 'fnscale' are supplied, 
  # optimization is carried out over the scaled function
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    parx = relist(flesh = parx*parscale, skeleton = skeleton)
    output = fn(parx, ...)/fnscale
    return(output)
  }
  
  # here, make methods explicit (one by one)
  output = .ahres(par=par, fn=fn1, gr=NULL, lower=lower, upper=upper, control=control, hessian=hessian)
  
  # reshaping full parameters, un-scaling par
  paropt = guess
  paropt[isActive] = output$par
  paropt = paropt*parscale
  output$par = paropt
  # value need to be unscaled
  output$value = output$value*fnscale
  
  output = c(output, active=list(par=isActive, flag=activeFlag))
  class(output) = c("ahres.results", class(output))
  
  return(output) 
  
}


# calibration_setup -------------------------------------------------------

#' Read and validate a calibration setup table
#'
#' Read a calibration settings table from a text file and coerce it to the
#' structure expected by [calibration_data()] and [calibration_objFn()]. The
#' setup table defines which model outputs are compared with observations, which
#' objective function is used for each output, whether a row contributes to the
#' calibration, and how each component is weighted in the final objective.
#'
#' The function is a validated wrapper around [utils::read.csv()]. It checks
#' required column names, coerces standard column types, fills a few optional
#' columns from `control`, and derives `weight` from `cv` when needed.
#'
#' @param file Path to the calibration setup file, typically a `.csv` file.
#'
#' @param control Named list of defaults used when optional columns are missing
#'   from the setup table. Currently recognised entries are:
#'   \describe{
#'     \item{`col_skip`}{Default number of leading columns to discard when a
#'     text file is read as a full table by [calibration_data()]. Defaults to
#'     `2`.}
#'     \item{`nrows`}{Default number of rows to read from each data file.
#'     Defaults to `NA`, which is interpreted downstream as “read all rows”.}
#'   }
#'
#' @param ... Additional arguments passed to [utils::read.csv()] when reading
#'   `file`.
#'
#' @return A `data.frame` describing the calibration targets. The returned table
#'   is suitable for direct use in [calibration_data()] and
#'   [calibration_objFn()]. Standard columns are coerced to the expected types
#'   and missing optional columns are added where needed.
#'
#' @details
#' The setup table may contain the following columns:
#'
#' \describe{
#'   \item{`variable`}{Name of the calibration component. This name must match
#'   the corresponding element returned by the model function supplied to
#'   [calibration_objFn()].}
#'
#'   \item{`type`}{Name of the objective function used for that component.
#'   In current workflows this is usually the name of a function accepted by
#'   [objFn()], such as `"norm2"`, `"lnorm2"`, `"pois"`, `"multinom"`,
#'   `"normp"`, or `"penalty"`.}
#'
#'   \item{`calibrate`}{Logical flag indicating whether the component is active
#'   in the calibration. Only rows with `calibrate = TRUE` are included in the
#'   objective function returned by [calibration_objFn()].}
#'
#'   \item{`weight`}{Positive weight used when partial objective values are
#'   aggregated. If `weight` is missing, or contains `NA`, it is derived from
#'   `cv` as `1 / (2 * cv^2)` for the corresponding rows.}
#'
#'   \item{`cv`}{Coefficient of variation used to derive `weight` when an
#'   explicit weight is not supplied. At least one of `weight` or `cv` must be
#'   available in the table.}
#'
#'   \item{`use_data`}{Logical flag indicating whether an observed dataset must
#'   be read for that row. This is typically `FALSE` for pure penalties such as
#'   `"normp"`, where the objective depends only on simulated values.}
#'
#'   \item{`file`}{Path to the observed data file for that component. Paths are
#'   interpreted relative to the `path` argument of [calibration_data()]. This
#'   column is still required when `use_data = FALSE`, although the file will not
#'   be read in that case.}
#'
#'   \item{`varid`}{Optional identifier used by [calibration_data()] when
#'   extracting a single column from a text file. If omitted, a full numeric
#'   table is read instead.}
#'
#'   \item{`col_skip`}{Optional number of leading columns to discard when a text
#'   file is read as a full table. Missing values are replaced by
#'   `control$col_skip`.}
#'
#'   \item{`nrows`}{Optional number of rows to read from each observed-data
#'   file. Missing values are replaced by `control$nrows`.}
#' }
#'
#' Required columns are `variable`, `type`, `calibrate`, `use_data`, and
#' `file`. In addition, at least one of `weight` or `cv` must be supplied.
#'
#' @section Weight handling:
#' Aggregation weights are stored in the returned table regardless of whether the
#' objective is later used in vector or scalar form.
#'
#' If the table does not contain a `weight` column, weights are computed from
#' `cv` as `1 / (2 * cv^2)`. If `weight` exists but contains missing values,
#' those entries are filled using the same formula. All final weights must be
#' strictly positive.
#'
#' @section Relationship with the rest of the workflow:
#' A typical file-based calibration workflow is:
#' \enumerate{
#'   \item read the setup table with `calibration_setup()`;
#'   \item load observed data with [calibration_data()];
#'   \item build the objective function with [calibration_objFn()];
#'   \item pass the resulting function to [calibrate()] or another optimiser.
#' }
#'
#' @examples
#' \dontrun{
#' demo = calibrar_demo(model = "PredatorPrey")
#' setup = calibration_setup(file = demo$setup)
#' setup
#' }
#'
#' @seealso [calibration_data()], [calibration_objFn()], [calibrate()], [objFn()]
#' @export
calibration_setup = function(file, control=list(), ...) {
  
  calibrationInfo = read.csv(file, stringsAsFactors=FALSE, ...)
  
  fullNames = c("variable", "type", "calibrate", "cv", "weight", "use_data", "file", "varid", "col_skip", "nrows")
  # mandatory columns
  minNames = c("variable", "type", "calibrate", "use_data", "file")
  doesNotMatch = !(names(calibrationInfo) %in% fullNames)
  dnm = names(calibrationInfo)[doesNotMatch]
  
  isMissing = !(minNames %in% names(calibrationInfo))
  im = minNames[isMissing]
  
  sdnm = if(length(dnm)>1) " columns do " else " column does "
  sim  = if(length(im)>1) " variables are " else " variable is "
  
  msg1 = sprintf("Error in %s (%s %s not match).", file, paste(sapply(dnm, sQuote), collapse=", "), sdnm)
  msg2 = sprintf("Error in %s (%s %s missing).", file, paste(sapply(im, sQuote), collapse=", "), sim)
  
  if(any(doesNotMatch)) stop(msg1)
  if(any(isMissing)) stop(msg2)
  
  if(is.null(calibrationInfo$weight) & is.null(calibrationInfo$cv))
    stop("Either 'cv' or 'weight' columns must be specified.")
  
  # parsing correct data types
  calibrationInfo$variable  = as.character(calibrationInfo$variable)
  calibrationInfo$type      = as.character(calibrationInfo$type)
  calibrationInfo$calibrate = as.logical(calibrationInfo$calibrate)
  if(!is.null(calibrationInfo$weight)) calibrationInfo$weight    = as.numeric(calibrationInfo$weight)
  if(!is.null(calibrationInfo$cv)) calibrationInfo$cv    = as.numeric(calibrationInfo$cv)
  calibrationInfo$use_data  = as.logical(calibrationInfo$use_data)
  calibrationInfo$file      = as.character(calibrationInfo$file)
  calibrationInfo$varid     = as.character(calibrationInfo$varid)
  
  if(is.null(control$col_skip)) control$col_skip = 2
  if(is.null(control$nrows)) control$nrows = NA

  # optional columns 
  
  if(is.null(calibrationInfo$weight))
    calibrationInfo$weight = 1/(2*calibrationInfo$cv^2)
  
  if(any(is.na(calibrationInfo$weight))) {
    if(is.null(calibrationInfo$cv)) stop("NAs found in 'weight' column.")
    isna = is.na(calibrationInfo$weight)
    calibrationInfo$weight[isna] = 1/(2*calibrationInfo$cv[isna]^2)
  }
  
  if(any(is.na(calibrationInfo$weight))) stop("NAs found in 'weight' column.")
  if(any(calibrationInfo$weight <= 0)) stop("All weights must be positive.")
  
  if(is.null(calibrationInfo$varid))
    calibrationInfo$varid = NA
  
  if(is.null(calibrationInfo$col_skip))
    calibrationInfo$col_skip = control$col_skip
  
  if(is.null(calibrationInfo$nrows))
    calibrationInfo$nrows = control$nrows
  
  return(calibrationInfo)
}

# calibration_data --------------------------------------------------------

#' Read observed data defined in a calibration setup table
#'
#' Build the `observed` object used by [calibration_objFn()] from a setup table
#' created with [calibration_setup()]. Each active row in the setup table points
#' to a file containing observed data, and this function reads those files into a
#' named list whose elements correspond to the calibration variables.
#'
#' Rows with `use_data = FALSE` are kept in the output as `NA`. This makes it
#' possible to include penalty-only components in the same calibration workflow:
#' such rows still participate in [calibration_objFn()] provided the model
#' returns a simulated component with the same name.
#'
#' @param setup A calibration setup table, usually created with
#'   [calibration_setup()].
#'
#' @param path Base directory used to resolve relative paths from `setup$file`.
#'   Each file is looked up with [base::file.path()] using `path` and the value
#'   stored in the corresponding `file` row.
#'
#' @param file Optional path to save the resulting `observed` object with
#'   [base::saveRDS()]. This is useful when observed data are expensive to read
#'   repeatedly.
#'
#' @param verbose Logical. If `TRUE`, print progress messages describing which
#'   variables were loaded.
#'
#' @param ... Additional arguments passed to [utils::read.csv()] when reading
#'   text files.
#'
#' @return A named list with one element per row of `setup`. Element names are
#'   taken from `setup$variable`.
#'
#'   Depending on the corresponding setup row, each element is typically one of:
#'   \itemize{
#'     \item a numeric vector, when `varid` selects a single column from a text
#'     file;
#'     \item a numeric matrix, when a full table is read from a text file;
#'     \item an arbitrary R object read from an `.rds` file;
#'     \item `NA`, when `use_data = FALSE` or `calibrate = FALSE`.
#'   }
#'
#' @details
#' Only rows satisfying both `setup$calibrate` and `setup$use_data` are actually
#' read from disk. All other rows are represented as `NA` in the returned list.
#'
#' File reading currently depends on file type:
#'
#' \describe{
#'   \item{text files (`.csv`, `.txt`)}{Read with [utils::read.csv()]. If
#'   `varid` is supplied, that column is extracted and returned as a numeric
#'   vector. Otherwise the full table is converted to a numeric matrix after
#'   dropping the first `col_skip` columns, when applicable.}
#'
#'   \item{R objects (`.rds`)}{Read with [base::readRDS()] and returned as-is.}
#'
#'   \item{netCDF (`.nc`)}{Not currently supported by the implementation.}
#' }
#'
#' The interpretation of `varid`, `col_skip`, and `nrows` is controlled by the
#' corresponding columns in `setup`, after defaults have been filled by
#' [calibration_setup()].
#'
#' @section Expected alignment with model outputs:
#' The names of the returned list are used later by [calibration_objFn()] to
#' match observed and simulated quantities. For each active calibration row,
#' `setup$variable` must therefore match the name of the corresponding element in
#' the list returned by the model function.
#'
#' @section Saving and reusing observed data:
#' When `file` is supplied, the created list is also saved as an `.rds` file.
#' This can be helpful in iterative calibration workflows where the setup table
#' is stable but reading the original data files repeatedly is inconvenient.
#'
#' @examples
#' \dontrun{
#' demo = calibrar_demo(model = "PredatorPrey")
#' setup = calibration_setup(file = demo$setup)
#' observed = calibration_data(setup = setup, path = demo$path)
#'
#' names(observed)
#' str(observed)
#' }
#'
#' @seealso [calibration_setup()], [calibration_objFn()], [objFn()]
#' @export
calibration_data = function(setup, path=".", file=NULL, verbose=TRUE, ...) {
  
  observed  = list()
  useData = as.logical(setup$use_data)
  isActive = as.logical(setup$calibrate)
  
  message("Creating observed data list for calibration...","\n")
  
  readData = useData & isActive
  
  for(iVar in seq_len(nrow(setup))) {
    
    ifile = file.path(path, setup$file[iVar])
    varid   = setup$varid[iVar]
    col_skip = setup$col_skip[iVar]
    nrows = setup$nrows[iVar]
    if(readData[iVar]) {
      observed[[iVar]] = .read_data(file=ifile, col_skip=col_skip, varid=varid, nrows=nrows, ...)
    } else {
      observed[[iVar]] = NA
    }
      
  } # end iVar loop
  
  msg0 = "Loaded observed data for variable: %s.\n"
  msg1 = "Loaded observed data for variables: %s.\n"
  
  msg = if(sum(readData, na.rm=TRUE)==1) msg0 else msg1 
  
  if(isTRUE(verbose)) message(sprintf(msg, paste(sQuote(setup$variable[readData]), collapse=", ")))
  
  names(observed) = as.character(setup$variable)
  
  if(!is.null(file)) saveRDS(observed, file=file)
  
  return(observed)
  
}

# calibration_objFn -------------------------------------------------------

#' Create a calibration objective function from a model, setup table, and data
#'
#' Build an objective function that runs a model, compares its outputs with the
#' observed data prepared by [calibration_data()], and returns either a vector of
#' partial objective values or a single aggregated objective value.
#'
#' This is the main bridge between a model-specific simulation function and the
#' optimisation tools in `calibrar`. The resulting function can be passed
#' directly to [calibrate()] and related optimisers.
#'
#' @param model Function that runs the model for a candidate parameter set. Its
#'   first argument must be `par`. The function must return a named list of
#'   simulated outputs whose names match `setup$variable` for all active
#'   calibration rows.
#'
#' @param setup A calibration setup table, usually created with
#'   [calibration_setup()].
#'
#' @param observed Observed data as returned by [calibration_data()]. This is
#'   typically a named list whose elements correspond to `setup$variable`.
#'
#' @param aggFn Function used to aggregate partial objective values when
#'   `aggregate = TRUE`. The default is an internal weighted sum. The function
#'   must accept arguments `x` (partial objective values) and `w` (weights).
#'
#' @param aggregate Logical. If `FALSE` (the default), the returned function
#'   produces one objective value per active calibration row. If `TRUE`, those
#'   partial values are combined into a single scalar using `aggFn` and the
#'   weights stored in `setup`.
#'
#' @param ... Additional arguments passed to `model` each time the returned
#'   objective function is evaluated.
#'
#' @return A function of `par` suitable for optimisation.
#'
#'   When evaluated, the returned function:
#'   \enumerate{
#'     \item runs `model(par, ...)`;
#'     \item extracts the simulated outputs;
#'     \item compares them with `observed` according to `setup$type`;
#'     \item returns either a named numeric vector of partial objective values
#'     or a single aggregated numeric value.
#'   }
#'
#'   The returned function has class `"objFn"` and carries metadata attributes
#'   used internally by [calibrate()], including the number of objective
#'   components, their weights, and the names of the active variables.
#'
#' @details
#' For each row where `setup$calibrate` is `TRUE`, the function looks for a
#' variable named `setup$variable` in both `observed` and the list returned by
#' `model`. The comparison itself is delegated to the objective function named in
#' `setup$type`, using the same convention as [objFn()].
#'
#' This means that the structure of the objects returned by `model` must already
#' be compatible with the corresponding objective function. For example:
#'
#' \itemize{
#'   \item pointwise criteria such as `"norm2"`, `"lnorm2"`, or `"pois"`
#'   usually expect vectors, matrices, or arrays with shapes matching the
#'   observed data;
#'   \item `"multinom"` expects matrix inputs;
#'   \item penalty terms such as `"normp"` may ignore `observed`, but still
#'   require the model output to contain the corresponding named element.
#' }
#'
#' If `model` returns `NULL`, the generated objective function returns `NULL`
#' after emitting a diagnostic message.
#'
#' @section Aggregated versus vector objectives:
#' With `aggregate = FALSE`, the returned objective is multi-component: it
#' produces one value per active calibration row. This can be useful for methods
#' that exploit vector-valued objectives.
#'
#' With `aggregate = TRUE`, the returned objective is scalar and can be used with
#' optimisers that require a single value. Aggregation is performed with
#' `aggFn(x = partial_values, w = setup$weight[setup$calibrate])`.
#'
#' @section Expected model interface:
#' In typical use, the model function has the form:
#'
#' `function(par, ...) -> named list`
#'
#' where `par` may be a numeric vector or a list, depending on how the
#' optimisation is configured. The returned list should expose one element per
#' calibration target, with names that match `setup$variable` exactly.
#'
#' @section Workflow:
#' The usual sequence is:
#' \enumerate{
#'   \item define a model-running function;
#'   \item read the setup table with [calibration_setup()];
#'   \item read observed data with [calibration_data()];
#'   \item create the objective function with `calibration_objFn()`;
#'   \item call [calibrate()] with the resulting function.
#' }
#'
#' @examples
#' \dontrun{
#' demo = calibrar_demo(model = "PredatorPrey")
#' setup = calibration_setup(file = demo$setup)
#' observed = calibration_data(setup = setup, path = demo$path)
#'
#' run_model = calibrar:::.PredatorPreyModel
#'
#' obj_vec = calibration_objFn(
#'   model = run_model,
#'   setup = setup,
#'   observed = observed,
#'   T = demo$T,
#'   aggregate = FALSE
#' )
#'
#' obj_scalar = calibration_objFn(
#'   model = run_model,
#'   setup = setup,
#'   observed = observed,
#'   T = demo$T,
#'   aggregate = TRUE
#' )
#'
#' obj_vec(demo$guess)
#' obj_scalar(demo$guess)
#' }
#'
#' @seealso [calibration_setup()], [calibration_data()], [calibrate()], [objFn()]
#' @export
calibration_objFn = function(model, setup, observed, aggFn=NULL, aggregate=FALSE, ...) {

  fn_name = deparse(substitute(model))
  
  fn   = match.fun(model)
  
  if(is.null(aggFn)) aggFn = .weighted.sum
  
  aggFn = match.fun(aggFn)
  
  force(observed)
  force(setup)
  force(aggregate)
  
  weights = setup$weight[setup$calibrate]

  msg = sprintf("The '%s' function returned NULL, please check.", fn_name)
  
  # check for names in observed and simulated
  fn1  = function(par) {
    aggFn = match.fun(aggFn)
    simulated = fn(par, ...)
    # in case the model produce a NULL, keep moving.
    if(is.null(simulated)) {
      message(msg)
      return(NULL)
    } 
    output = .calculateObjetiveValue(obs=observed, sim=simulated, info=setup)
    if(isTRUE(aggregate)) output = aggFn(x=output, w=weights)
    return(output)
  }
  
  class(fn1) = c(class(fn1), "objFn")

  fnx = function(par) {
    simulated = fn(par, ...)
    return(simulated)
  }
  
  attr(fn1, "nvar") = if(aggregate) 1 else sum(setup$calibrate)
  attr(fn1, "weights") = weights
  attr(fn1, "variables") = setup$variable[setup$calibrate]
  attr(fn1, "aggregate") = aggregate
  attr(fn1, "fn") = fnx
  attr(fn1, "aggFn") = if(aggregate) NULL else aggFn
  return(fn1) 
  
}

multiMethods = "AHR-ES" # list of methods supporting multi-objective
