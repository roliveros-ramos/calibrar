# getObservedData ---------------------------------------------------------

#' @title Get observed data for the calibration of a model 
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
  
  .Deprecated("calibration_data", package="calibrar")
  
  observed  = list()
  variables = info$variable
  
  useData = as.logical(info$useData)
  
  message("Creating observed data list for calibration...","\n")
  
  for(iVar in seq_len(nrow(info))) {
    
    message(paste0("Variable: ", variables[iVar], "\n"))
    varPath         = file.path(path, data.folder, paste0(variables[iVar],".csv"))
    observed[[iVar]] = if(useData[iVar]) .read.csv3(varPath, ...) else NA
    
  }
  
  names(observed) = variables
  
  return(observed)
  
}

# getCalibrationInfo ------------------------------------------------------

#' Get information to run a calibration using the \code{calibrar} package.
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
getCalibrationInfo = function(path, file="calibrationInfo.csv", stringsAsFactors=FALSE, ...) {
  
  .Deprecated("calibration_setup", package="calibrar")
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

#' Create an objective function to be used with optimization routines
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
  
  .Deprecated("calibration_objFn", package="calibrar")
  
  fn   = match.fun(runModel)
  aggFn = match.fun(aggFn)
  
  force(observed)
  force(info)
  force(aggregate)
  
  weights = info$weights[info$calibrate]
  
  # check for names in observed and simulated
  fn1  = function(par) {
    aggFn = match.fun(aggFn)
    simulated = fn(par, ...)
    # apply fitness to all outputs
    output = .calculateObjetiveValue(obs=observed, sim=simulated, info=info)
    if(isTRUE(aggregate)) output = aggFn(x=output, w=weights)
    return(output)
  }
  
  class(fn1) = c(class(fn1), "objFn")
  
  fnx = function(par) {
    simulated = fn(par, ...)
    return(simulated)
  }
  
  attr(fn1, "nvar") = sum(info$calibrate)
  attr(fn1, "weights") = weights
  attr(fn1, "variables") = info$variables[info$calibrate]
  attr(fn1, "aggregate") = aggregate
  attr(fn1, "fn") = fnx
  return(fn1) 
  
}

#' @rdname calibrar_demo
#' @export 
calibrarDemo = function(path=NULL, model=NULL,  ...) {
  
  .Deprecated("calibrar_demo", package="calibrar")
  return(calibrar_demo(path=path, model=model,  ...))                
  
}