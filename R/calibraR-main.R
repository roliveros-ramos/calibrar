# calibrar package: Automated Calibration for Complex (Ecological) Models --------

#' @title Automated Calibration for Complex (Ecological) Models
#' @description Automated Calibration for Complex (Ecological) Models
#' @name calibrar-package
#' @aliases calibrar-package calibrar
#' @docType package
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
#' demo = calibrarDemo(model="PoissonMixedModel", L=5, T=100) 
#' # get calibration information
#' calibrationInfo = getCalibrationInfo(path=demo$path)
#' # get observed data
#' observed = getObservedData(info=calibrationInfo, path=demo$path)
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
#' obj  = createObjectiveFunction(runModel=runModel, info=calibrationInfo, 
#'                                observed=observed, forcing=forcing)
#' cat("Starting calibration...\n")
#' control = list(weights=calibrationInfo$weights, maxit=3.6e5) # control parameters
#' cat("Running optimization algorithms\n", "\t", date(), "\n")
#' cat("Running optim AHR-ES\n")
#' ahr = calibrate(par=demo$guess, fn=obj, lower=demo$lower, upper=demo$upper, control=control)
#' summary(ahr)
#' } 
NULL

# calibrate ---------------------------------------------------------------

#' @title Sequential parameter estimation for the calibration of models
#' @description This function performs the optimization of a function, possibly 
#' in sequential phases of increasing complexity, and it is designed for the 
#' calibration of a model, by minimizing the error function \code{fn} associated to it.  
#' @param par A numeric vector. The length of the par argument defines the 
#' number of parameters to be estimated (i.e. the dimension of the problem).
#' @param fn The function to be minimized.
#' @param gr the gradient of \code{fn}. Ignored, added for portability with
#' other optimization functions.
#' @param \dots Additional parameters to be passed to \code{fn}.
#' @param method The optimization method to be used. The default method
#' is the AHR-ES (Adaptative Hierarchical Recombination Evolutionary Strategy, 
#' Oliveros-Ramos & Shin, 2016). All the methods from stats::optim,
#' optimx::optimx and cmaes::cma_es are available.
#' @param lower Lower threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{-Inf}. By default \code{-Inf} is used (unconstrained).
#' @param upper Upper threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{Inf}. By default \code{Inf} is used (unconstrained). 
#' @param control Parameter for the control of the algorithm itself, see details.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Currently not implemented. 
#' @param phases An optional vector of the same length as \code{par}, 
#' indicating the phase at which each parameter becomes active. If omitted, 
#' default value is 1 for all parameters, performing a single optimization.
#' @param replicates The number of replicates for the evaluation of \code{fn}.
#' The default value is 1. A value greater than 1 is only useful for stochastic
#' functions.
#' @details In the control list, \code{aggFn} is a function to aggregate \code{fn} to 
#' a scalar value if the returned value is a vector. Some optimization algorithm can 
#' exploite the additional information provided by a vectorial output from \code{fn}.
#' @author Ricardo Oliveros-Ramos
#' @examples
#' calibrate(par=rep(NA, 5), fn=SphereN)
#' \dontrun{
#' calibrate(par=rep(NA, 5), fn=SphereN, replicates=3)
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5)
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5, phases=c(1,1,1,2,3))
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=c(1,1,4), lower=-5, upper=5, phases=c(1,1,1,2,3))
#' }
#' @export
calibrate = function(par, fn, gr, ..., method, lower, upper, control, 
                     hessian, phases, replicates) {
  UseMethod("calibrate")
}


#' @export
calibrate.default = function(par, fn, gr = NULL, ..., method = "AHR-ES",
                     lower = NULL, upper = NULL, control = list(), 
                     hessian = FALSE, phases = NULL, replicates=1) {

  if(method=="AHR-ES") method = "default"
  
  # check function and method
  multiMethods = "default" # list of methods supporting multi-objective
  
  if(inherits(fn, "objFn") & !(method %in% multiMethods)) {
    agg = attr(fn, "aggregate")
    if(is.null(agg)) warning("Update your objective function to the last version of the package.")
    if(!isTRUE(agg)) 
      stop(sprintf("Method '%s' does not support multi-objective optimization, use aggregate=TRUE in 'createObjectiveFunction'.",
                   method))
  }
  
  # check for a restart file
  restart = .restartCalibration(control, type="results")
  
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
  
  output = if(isTRUE(restart)) .getResults(control=control) else list(phase=1)
  
  conv = .checkConvergence(control, nphases)
  
  # start the sequential parameter estimation
  for(phase in seq(from=output$phase, to=nphases)) {
    
    if(output$phase > nphases) break
  
    control$maxgen      = conv$maxgen[phase]
    control$maxiter     = conv$maxiter[phase]
    control$convergence = conv$convergence[phase]

    active = (phases <= phase) # NAs are corrected in .calibrar 
    # call optimizers handler .calibrar
    
    tm1 = Sys.time()
    temp = .calibrar(par=par, fn=fn, gr = gr, ..., method = method, 
                   lower = lower, upper = upper, control = control, 
                   hessian = hessian, active=active, skeleton=skeleton,
                   replicates=replicates[phase])
    tm2 = Sys.time()
    
    output$phases[[phase]] = temp # trim?
    output$phase = phase + 1
    
    .createOutputFile(output, control) 
    
    par = temp$par #update parameter guess
    control = .updateControl(control=control, opt=temp, method=method)  # update CVs? 
    
    msg = paste(c(sprintf("\nPhase %d finished in %s (%d of %d parameters active)",
                        phase, format_difftime(tm1, tm2), sum(active, na.rm=TRUE), npar),
                sprintf("Function value: %g", temp$value),
                paste(c("Parameter values:",sprintf("%0.3g", par[which(active)])), collapse=" "), 
                "\n"), collapse="\n")
    message(msg)
  }
  
   isActive = (phases>0) & !is.na(phases)
   paropt = output$phases[[nphases]]$par # save parameters of last phase
  
#   newNames = rep("*", npar)
#   newNames[isActive] = ""
#   
#   names(paropt) = paste0(names(paropt), newNames)

  paropt = relist(paropt, skeleton)
  class(paropt) = setdiff(class(paropt), "relistable")
  
  final = list(par=paropt, value=output$phases[[nphases]]$value, 
               counts=output$phases[[nphases]]$counts, 
               partial=output$phases[[nphases]]$partial, 
               active=isActive, fn=fn)
  
  output = c(final, output)
  class(output) = c("calibrar.results")
  .createOutputFile(output, control) 
  
  return(output)
  
}

# optimES -----------------------------------------------------------------

#' @title Optimization using Evolutionary Strategies
#' @description This function performs the optimization of a function using 
#' evolutionary strategies, by default the AHR-ES (Oliveros & Shin, 2015). 
#' @param par A numeric vector. The length of the par argument defines the 
#' number of parameters to be estimated (i.e. the dimension of the problem).
#' @param fn The function to be minimized.
#' @param gr the gradient of \code{fn}. Ignored, added for portability with
#' other optimization functions.
#' @param \dots Additional parameters to be passed to \code{fn}.
#' @param lower Lower threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{-Inf}. By default \code{-Inf} is used (unconstrained).
#' @param upper Upper threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{Inf}. By default \code{Inf} is used (unconstrained). 
#' @param active A boolean vector of the same length of \code{par}. If \code{TRUE}, the
#' parameter is optimized, if \code{FALSE} the parameter is fixed to the value specified
#' in \code{par}.
#' @param control Parameter for the control of the algorithm itself, see details.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Currently not implemented. 
#' @param method The optimization method to be used. Currently, the only implemented
#' is the 'default' method, corresponding to the AHR-ES (Oliveros & Shin, 2015).
#' @author Ricardo Oliveros-Ramos
#' @examples
#' optimES(par=rep(1, 5), fn=SphereN)
#' @export
optimES = function (par, fn, gr = NULL, ..., lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, method = "default") {
  
  
  npar = length(par)

  if(is.null(names(par))) names(par) = .printSeq(npar, preffix="par")
  
  active = .checkActive(active=active, npar=npar)
  bounds = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess  = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  par    = guess[isActive]
  lower  = bounds$lower[isActive]
  upper  = bounds$upper[isActive]
  
  npar = length(par)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  
  control = .checkControl(control=control, method=method, par=par, fn=fn, active=active, ...)
  
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    fn(parx, ...)/control$fnscale
  }
  
  output = .optimES(par=par, fn=fn1, lower=lower, upper=upper, control=control, isActive = isActive)
  
  paropt = guess
  paropt[isActive] = output$ppar 
  
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  output = list(par=paropt, output, active=list(par=isActive, flag=activeFlag))
  
  class(output) = c("optimES.result", class(output))
  
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
  
  fullNames = c("variable", "type", "calibrate", "weight", "use_data", "file", "varid")  
  doesNotMatch = !(names(calibrationInfo) %in% fullNames)
  dnm = names(calibrationInfo)[doesNotMatch]
  
  isMissing = !(fullNames %in% names(calibrationInfo))
  im = fullNames[isMissing]
  
  sdnm = if(length(dnm)>1) " columns do " else " column does "
  sim  = if(length(im)>1) " variables are " else " variable is "
  
  msg1 = sprintf("Error in %s (%s %s not match).", file, paste(sapply(dnm, sQuote), collapse=", "), sdnm)
  msg2 = sprintf("Error in %s (%s %s not match).", file, paste(sapply(im, sQuote), collapse=", "), sim)
  
  if(any(doesNotMatch)) stop(msg1)
  if(any(isMissing)) stop(msg2)
  
  # cating correct data types
  calibrationInfo$variable  = as.character(calibrationInfo$variable)
  calibrationInfo$type      = as.character(calibrationInfo$type)
  calibrationInfo$calibrate = as.logical(calibrationInfo$calibrate)
  calibrationInfo$weight    = as.numeric(calibrationInfo$weight)
  calibrationInfo$use_data  = as.logical(calibrationInfo$use_data)
  calibrationInfo$file      = as.character(calibrationInfo$file)
  calibrationInfo$varid     = as.character(calibrationInfo$varid)
  
  if(is.null(control$col_skip)) control$col_skip = 2
    
  if(is.null(calibrationInfo$col_skip))
    calibrationInfo$col_skip = control$col_skip
  
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
#' @param \dots Additional arguments to \code{read.csv} function 
#' to read the data files.
#' @return A list with the observed data needed for a calibration, to be used 
#' in combination with the \code{\link{createObjectiveFunction}}.
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{createObjectiveFunction}}, \code{\link{getCalibrationInfo}}.
#' @export
calibration_data = function(setup, path=".", verbose=TRUE, ...) {
  
  observed  = list()
  useData = as.logical(setup$use_data)
  
  message("Creating observed data list for calibration...","\n")
  
  for(iVar in seq_len(nrow(setup))) {
    
    if(isTRUE(verbose)) message(sprintf("Variable: %s", setup$variable[iVar]))
    ifile = file.path(path, setup$file[iVar])
    varid   = setup$varid[iVar]
    col_skip = setup$col_skip[iVar]
    if(useData[iVar]) {
      observed[[iVar]] = .read.csv4(file=ifile, col_skip=col_skip, varid=varid, ...)
    } else {
      observed[[iVar]] = NA
    }
      
  }
  
  names(observed) = as.character(setup$variable)
  
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
calibration_objFn = function(model, setup, observed, aggFn=NULL, 
                                   aggregate=FALSE, ...) {

  fn   = match.fun(model)
  
  if(is.null(aggFn)) aggFn = .weighted.sum
  
  aggFn = match.fun(aggFn)
  
  force(observed)
  force(setup)
  force(aggregate)
  
  weights = setup$weight[setup$calibrate]

  # check for names in observed and simulated
  fn1  = function(par) {
    aggFn = match.fun(aggFn)
    simulated = fn(par, ...)
    # apply objFn to all outputs
    output = .calculateObjetiveValue(obs=observed, sim=simulated, info=setup)
    if(isTRUE(aggregate)) output = aggFn(x=output, w=weights)
    return(output)
  }
  
  class(fn1) = c(class(fn1), "objFn")

  fnx = function(par) {
    simulated = fn(par, ...)
    return(simulated)
  }
  
  attr(fn1, "nvar") = sum(setup$calibrate)
  attr(fn1, "weights") = weights
  attr(fn1, "variables") = setup$variable[setup$calibrate]
  attr(fn1, "aggregate") = aggregate
  attr(fn1, "fn") = fnx
  return(fn1) 
  
}
