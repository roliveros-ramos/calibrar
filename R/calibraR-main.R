# calibrar package: Calibration of Ecological models using Evoluti --------

#' Calibration of Ecological Models using Evolutionary Algorithms
#' 
#' Calibration of Ecological Models using Evolutionary Algorithms
#' 
#' \tabular{ll}{ Package: \tab calibrar\cr Type: \tab Package\cr Version: \tab
#' 0.1\cr Date: \tab 2014-09-14\cr License: \tab GPL-2\cr } calibrate()
#' optimES()
#' 
#' @name calibrar-package
#' @aliases calibrar-package calibrar
#' @docType package
#' @author Ricardo Oliveros-Ramos Maintainer: Ricardo Oliveros-Ramos
#' <ricardo.oliveros@@gmail.com>
#' @references calibrar: an R package for the calibration of ecological models (Oliveros-Ramos and Shin 2014)
#' @keywords calibration
#' @examples
#' 
#' calibrate(par=rep(NA, 5), fn=SphereN)
#' calibrate(par=rep(NA, 5), fn=SphereN, replicates=3)
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5)
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5, phases=c(1,1,1,2,3))
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=c(1,1,4), lower=-5, upper=5, phases=c(1,1,1,2,3))
#' 
NULL

# optimES -----------------------------------------------------------------

#' @title Optimization using Evolutionary Strategies
#' @description This function performs the optimization of a function using 
#' evolutionary strategies, by default the AHR-ES (Oliveros & Shin, 2014). 
#' @param par A numeric vector. The length of the par argument defines the 
#' number of parameters to be estimated (i.e. the dimension of the problem).
#' @param fn The function to be minimized.
#' @param gr the gradient of \code{fn}. Ignored, added for portability with
#' other optimization functions.
#' @param \dots Additional parametrs to be passed to \code{fn}.
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
#' is the 'default' method, corresponding to the AHR-ES (Oliveros & Shin, 2014).
#' @param restart If restart is an object of class \code{optimES.restart}, it continue the
#' optimization using the information from that object.
#' @author Ricardo Oliveros-Ramos
#' @examples
#' calibrate(par=rep(1, 5), fn=SphereN)
#' @export
optimES = function (par, fn, gr = NULL, ..., lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, method = "default", 
                    restart=NULL) {
  
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
  
  control = .checkControl(control=control, method=method, par=par, fn=fn, active=active)
  
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    fn(parx, ...)/control$fnscale
  }
  
  trace = "elmo"
  
  if(control$REPORT>0 & control$trace>0) {
    
    trace = list()
    trace$control = control
    trace$par = matrix(NA, nrow=control$maxgen, ncol=length(isActive))
    trace$value = rep(NA, control$maxgen)
    trace$best  = rep(NA, control$maxgen)
    
    if(control$trace>1) {
      trace$sd = matrix(NA, nrow=control$maxgen, ncol=length(isActive))   
      trace$step = rep(NA, control$maxgen)     
    }

    if(control$trace>2) trace$opt = vector("list", control$maxgen)
    
  } 

  
  # opt = get restart, a method for a file (character) or a restart class
  
  opt = if(!is.null(restart)) # TO_DO
    .getRestart(restart=restart) else 
      .newOpt(par=par, lower=lower, upper=upper, control=control) # here par=NA?
  
  while(isTRUE(.continueEvolution(opt, control))) {
    
    opt$gen  = opt$gen + 1
    opt$ages = opt$ages + 1
    
    # create a new population
    
    if(all(opt$SIGMA==0)) break
    opt$pop = .createPopulation(opt)
    
    # evaluate the function in the population: evaluate fn, aggregate fitness
    
    opt$fitness = .calculateFitness(opt, fn=fn1)
    
    # select best 'individuals'
    
    opt$selected = .selection(opt)
        
    # create the new parents: MU and SD
    
    opt = .calculateOptimalWeights(opt)
    opt = .updatePopulation(opt)
    
    # save status of the population (for restart)
    
    # save detailed outputs
    if(control$REPORT>0 & control$trace>0) {
        
        trace$par[opt$gen, ] = opt$MU
        trace$best[opt$gen]  = opt$selected$best$fit.global
        
        if(control$trace>1) {
          trace$sd[opt$gen, ]  = opt$SIGMA
          trace$step[opt$gen]  = opt$step       
        }

      if(opt$gen%%control$REPORT==0) {
        trace$value[opt$gen] = control$aggFn(fn1(opt$MU), control$weights)
        if(control$trace>2) trace$opt[[opt$gen]] = opt
      }
      
    }
    
   if(control$verbose & opt$gen%%control$REPORT==0) 
     .messageByGen(opt, trace)
    
  } # end generations loop
  
  value = control$aggFn(x=fn1(opt$MU), w=control$weights)
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)
  
  paropt = guess
  paropt[isActive] = opt$MU
  
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  newNames = rep("*", npar)
  newNames[isActive] = ""
  
  names(paropt) = paste0(names(paropt), newNames)
  
  output = list(par=paropt, value=value, counts=opt$counts, 
                trace=trace, partial=fn1(opt$MU), MU=opt$MU, 
                active=list(par=isActive, flag=activeFlag))
  
  class(output) = c("optimES.result", class(output))
  
  return(output)
  
}

# calibrate ---------------------------------------------------------------

#' @title Sequential calibration of models
#' @description This function performs the optimization of a function, possibly 
#' in sequential phases of increasing complexity, and it is designed for the 
#' calibration of a model, by minimizing the error function \code{fn} associated to it.  
#' @param par A numeric vector. The length of the par argument defines the 
#' number of parameters to be estimated (i.e. the dimension of the problem).
#' @param fn The function to be minimized.
#' @param \dots Additional parametrs to be passed to \code{fn}.
#' @param aggFn A function to aggregate \code{fn} to a scalar value if the
#' returned value is a vector. Some optimization algorithm can explote the
#' additional information provided by a vectorial output from \code{fn}.
#' @param phases An optional vector of the same length as \code{par}, 
#' indicating the phase at which each parameter becomes active. If omitted, 
#' default value is 1 for all parameters, performing a single optimization.
#' @param replicates The number of replicates for the evaluation of \code{fn}.
#' The default value is 1. A value greater than 1 is only useful for stochastic
#' functions.
#' @param lower Lower threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{-Inf}. By default \code{-Inf} is used (unconstrained).
#' @param upper Upper threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{Inf}. By default \code{Inf} is used (unconstrained). 
#' @param gr the gradient of \code{fn}. Ignored, added for portability with
#' other optimization functions.
#' @param control Parameter for the control of the algorithm itself, see details.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Currently not implemented. 
#' @param method The optimization method to be used. Currently, the only implemented
#' is the 'default' method, corresponding to the AHR-ES (Oliveros & Shin, 2014).
#' @param restart If restart is an object of class \code{calibrar.restart}, it continue the
#' calibration using the information from that object.
#' @author Ricardo Oliveros-Ramos
#' @examples
#' calibrate(par=rep(1, 5), fn=SphereN)
#' @export
calibrate = function(par, fn, ..., aggFn = NULL, phases = NULL, replicates=1, 
                     lower = -Inf, upper = Inf,  gr = NULL, control = list(), 
                     hessian = FALSE, method = "default", restart = NULL) {
  
  npar = length(par)
  
  fn = match.fun(fn)
  
  phases     = .checkPhases(phases=phases, npar=npar)
  bounds     = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess      = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)

  par     = guess
  lower   = bounds$lower
  upper   = bounds$upper
  nphases = max(phases, na.rm=TRUE)

  replicates = .checkReplicates(replicates, nphases) 
  
  output = list()
  
  for(phase in seq_len(nphases)) {
      
    # call optimEA
    active = (phases <= phase) # NAs are corrected in optimES 
    temp = optimES(par=par, fn=fn, gr = NULL, ..., method = method, 
                   lower = lower, upper = upper, active=active, 
                   control = control, hessian = hessian, restart=restart)
   
    output$phases[[phase]] = temp # trim?
    par[which(active)] = temp$MU
    control = .updateControl(control=control, opt=temp, method=method)  # update CVs? 

    cat(sprintf("\nPhase %d finished (%d of %d parameters active)\n",
                phase, sum(active, na.rm=TRUE), npar))
    cat(sprintf("Function value: %g \n", temp$value))
    print(temp$MU)
    cat("\n")
  }
  
  isActive = !is.na(phases) & (phases>=1)
  paropt = guess
  paropt[isActive] = temp$MU
  
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  newNames = rep("*", npar)
  newNames[isActive] = ""
  
  names(paropt) = paste0(names(paropt), newNames)
  
  final = list(par=paropt, value=temp$value, counts=temp$counts, 
               partial=temp$partial, active=isActive)
  
  output = c(final, output)
  class(output) = c("calibrar.result", class(output))
  
  return(output)
  
}

# getObservedData ---------------------------------------------------------

#' @title Get observed data for a model calibration 
#' 
#' @description Create a list with the observed data with the 
#' information provided by its main argument. 
#' 
#' @param info A data.frame with the information about the calibration, 
#' normally created with the \code{\link{getCalibrationInfo}} function. 
#' See details.
#' @param path Path to the directory to look up for the data.
#' @param data.folder folder in the path containing the data.
#' @param \dots Additional arguments to \code{read.csv} function 
#' to read the data files.
#' @return A list with the observed data needed for a calibration, to be used 
#' in combination with the \code{\link{createObjectiveFunction}}.
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{createObjectiveFunction}}, \code{\link{getCalibrationInfo}}.
#' @export
getObservedData = function(info, path, data.folder="data", ...) {
  
  observed  = list()
  variables = info$variable
  
  useData       = as.logical(info$useData)
  
  cat("Creating observed data list for calibration...","\n")
  
  for(var in 1:nrow(info)) {
    
    cat(paste0("Variable: ", variables[var], "\n"))
    var.path        = file.path(path, data.folder, paste0(variables[var],".csv"))
    datos           = if(useData[var]) .read.csv3(var.path, ...) else NA
    observed[[var]] = datos
    
  }
  
  names(observed) = variables
  
  return(observed)
  
}

# getCalibrationInfo ------------------------------------------------------

#' Get information for a calibration using the \code{calibrar} package.
#' 
#' A wrapper for \code{read.csv} checking column names and data types 
#' for the table with the calibration information.
#' 
#' @param path The path to look for the file.
#' @param file The file with the calibration information, see details.
#' @param stringsAsFactors To be passed to \code{read.csv}.
#' @param \dots Additional arguments to \code{read.csv} function.
#' @return A data.frame with the information for the calibration of a 
#' model, to be used with the \code{\link{createObjectiveFunction}} 
#' and \code{\link{getObservedData}}.
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{createObjectiveFunction}}, \code{\link{getObservedData}}.
#' @export
getCalibrationInfo = function(path, file="calibrationInfo.csv", 
                              stringsAsFactors=FALSE, ...) {
  
  caliPath = file.path(path, file)
  calibrationInfo = read.csv(caliPath, stringsAsFactors=FALSE, ...)
  
  fullNames = c("variable", "type", "calibrate", "weights", "useData")  
  doesNotMatch = !(names(calibrationInfo) %in% fullNames)
  dnm = names(calibrationInfo)[doesNotMatch]
  
  isMissing = !(fullNames %in% names(calibrationInfo))
  im = fullNames[isMissing]
  
  sdnm = if(length(dnm)>1) " columns do " else " column does "
  sim  = if(length(im)>1) " variables are " else " variable is "
  msg1 = paste0("Error in ", caliPath, " file (", paste(sapply(dnm, sQuote), collapse=", "), 
                sdnm, "not match).")
  msg2 = paste0("Error in ", caliPath, " file (", paste(sapply(im, sQuote), collapse=", "), 
                sim, "missing).")
  
  if(any(doesNotMatch)) stop(msg1)
  if(any(isMissing)) stop(msg2)
  
  # cating correct data types
  calibrationInfo$variable  = as.character(calibrationInfo$variable)
  calibrationInfo$type      = as.character(calibrationInfo$type)
  calibrationInfo$calibrate = as.logical(calibrationInfo$calibrate)
  calibrationInfo$weights   = as.numeric(calibrationInfo$weights)
  calibrationInfo$useData   = as.logical(calibrationInfo$useData)
  
  return(calibrationInfo)
}

# createObjectiveFunction -------------------------------------------------

#' Create and objective function to be used with optimization routines
#' 
#' Create a new function, to be used as the objective function in the 
#' calibration, given a function to run the model within R, observed data 
#' and information about the comparison with data.
#' 
#' @param runModel Function to run the model and produce a list of outputs.
#' @param info A data.frame with the information about the calibration, 
#' normally created with the \code{\link{getCalibrationInfo}} function. 
#' See details.
#' @param observed A list of the observed variables created with the 
#' function \code{\link{getObservedData}}
#' @param aggFn A function to aggregate \code{fn} to a scalar value if the
#' returned value is a vector. Some optimization algorithm can explote the
#' additional information provided by a vectorial output from \code{fn}
#' @param aggregate boolean, if TRUE, a scalar value is returned using the 
#' \code{aggFn}.
#' @param \dots More arguments passed to the \code{runModel} function.
#' @return A function, integrating the simulation of the model and the 
#' comparison with observed data. 
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{getObservedData}}, \code{\link{getCalibrationInfo}}.
#' @export
createObjectiveFunction = function(runModel, info, observed, aggFn=.weighted.sum, 
                                   aggregate=FALSE, ...) {

  fn   = match.fun(runModel)
  aggFn = match.fun(aggFn)
  
  force(observed)
  force(info)
  force(aggregate)
  
  # check for names in observed and simulated
  
  fn1  = function(par) {
    aggFn = match.fun(aggFn)
    simulated = fn(par, ...)
    # apply fitness to all outputs
    output = .calculateObjetiveValue(obs=observed, sim=simulated, info=info)
    if(isTRUE(aggregate)) output = aggFn(x=output, w=info$weights)
    return(output)
  }
  
  return(fn1) 
  
}
