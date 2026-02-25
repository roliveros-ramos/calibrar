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

#' General-purpose optimization using heuristic algorithms
#'
#' @param gr Function to compute the gradient of \code{fn}. Ignored by most methods, 
#' added for consistency with other optimization functions.
#' @examples 
#' optim2(par=rep(NA, 5), fn=sphereN)
#' @inheritParams optim2
#' @inherit optim2 return author
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

# AHR-ES ------------------------------------------------------------------

#' @title Adaptative Hierarchical Recombination Evolutionary Strategy (AHR-ES) for 
#' derivative-free and black-box optimization 
#' @description This function performs the optimization of a function using 
#' the Adaptative Hierarchical Recombination Evolutionary Strategy (AHR-ES, Oliveros & Shin, 2015). 
#'
#' @author Ricardo Oliveros-Ramos
#' @examples
#' \dontrun{ahres(par=rep(1, 5), fn=sphereN)}
#' @inheritParams calibrate
#' @inheritParams optim2
#' @inherit optim2 return
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
  class(output) = c("ahres.result", class(output))
  
  return(output) 
  
}


# calibration_setup -------------------------------------------------------

#' Get information to run a calibration using the \code{calibrar} package.
#' 
#' A wrapper for \code{read.csv} checking column names and data types 
#' for the table with the calibration information.
#' 
#' @param file The file with the calibration information, see details.
#' @param control Control arguments for generating the setup. See details.
#' @param \dots Additional arguments to \code{read.csv} function.
#' @return A data.frame with the information for the calibration of a 
#' model, to be used with the \code{\link{calibration_objFn}} 
#' and \code{\link{calibration_data}}.
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{calibration_objFn}}, \code{\link{calibration_data}}.
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

#' @title Get observed data for the calibration of a model 
#' 
#' @description Create a list with the observed data with the 
#' information provided by its main argument. 
#' 
#' @param setup A data.frame with the information about the calibration, 
#' normally created with the \code{\link{calibration_setup}} function. 
#' See details.
#' @param path Path to the directory to look up for the data. Paths in setup are considered relatives to this path.
#' @param file Optional file to save the created object (as an 'rds' file.)
#' @param verbose If TRUE, detailed messages of the process are printed.
#' @param \dots Additional arguments to \code{read.csv} function 
#' to read the data files.
#' @return A list with the observed data needed for a calibration, to be used 
#' in combination with the \code{\link{calibration_objFn}}.
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{calibration_objFn}}, \code{\link{calibration_setup}}.
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

#' Create an objective function to be used with optimization routines
#' 
#' Create a new function, to be used as the objective function in the 
#' calibration, given a function to run the model within R, observed data 
#' and information about the comparison with data.
#' 
#' @param model Function to run the model and produce a list of outputs.
#' @param setup A data.frame with the information about the calibration, 
#' normally created with the \code{\link{calibration_setup}} function. 
#' See details.
#' @param observed A list of the observed variables created with the 
#' function \code{\link{calibration_data}}
#' @param aggFn A function to aggregate \code{fn} to a scalar value if the
#' returned value is a vector. Some optimization algorithm can explote the
#' additional information provided by a vectorial output from \code{fn}
#' @param aggregate boolean, if TRUE, a scalar value is returned using the 
#' \code{aggFn}.
#' @param \dots More arguments passed to the \code{model} function.
#' @return A function, integrating the simulation of the model and the 
#' comparison with observed data. 
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{calibration_data}}, \code{\link{calibration_setup}}.
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
