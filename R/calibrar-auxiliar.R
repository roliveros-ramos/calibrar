
# calibrarDemo ------------------------------------------------------------

#' @title Demos for the calibrar package
#' @description Creates demo files able to be processed for a full calibration using
#' the calibrar package
#' 
#' @param path Path to create the demo files 
#' @param model Model to be used in the demo files, see details. 
#' @param \dots Additional parameters to be used in the construction of
#' the demo files.
#' @return A list with the following elements:
#' \item{path}{Path were the files were saved}
#' \item{par}{Real value of the parameters used in the demo} 
#' \item{constants}{Constants used in the demo} 
#' @author Ricardo Oliveros--Ramos
#' @references Oliveros-Ramos and Shin (2014)
#' @keywords demo calibration 
#' @examples
#' \dontrun{
#' require(calibrar)
#' path = NULL #' NULL to use the current directory
#' 
#' # create the demonstration files
#' demo = calibrarDemo(path=path, model="PoissonMixedModel", L=4, T=20) 
#' 
#' # set.seed(12345) #' updated to T=20 and L=40 for comparative purposes.
#' # Parameter information
#' parInfo = read.csv(file.path(demo$path, "parInfo.csv"), row.names=1)
#' 
#' # get calibration information
#' calibrationInfo = getCalibrationInfo(path=demo$path)
#' 
#' # get observed data
#' observed = getObservedData(info=calibrationInfo, path=demo$path)
#' 
#' # read forcings for the model
#' forcing = read.csv(file.path(demo$path, "master", "environment.csv"),
#'                    row.names=1)
#' 
#' # Defining 'runModel' function
#' runModel = function(par, forcing) {
#'   
#'   #' forcing is a matrix with the values of environmental variables
#'   T = nrow(forcing) #' get number of time steps
#'   L = ncol(forcing) #' get number of sites
#'   
#'   # create parameter list in a format readable by model
#'   parList = list()
#'   parList$alpha = par[1]
#'   parList$beta  = par[2]
#'   parList$gamma = par[2 + seq_len(T-1)]
#'   parList$sd = par[T+2]
#'   parList$mu_ini = par[T + 2 + seq_len(L)]
#'   
#'   # get the model
#'   model = calibrar:::.PoissonMixedModel
#'   # run the model
#'   output = model(par=parList, forcing=forcing)
#'   output = as.list(as.data.frame(output)) #' return a list with the results
#'   # names of the outputs matching observed data names
#'   names(output) = paste0("site_", sprintf(paste0("%0", ceiling(log10(L+1)), "d"), seq_len(L)))
#'   
#'   output = c(output, list(gammas=parList$gamma))
#'   return(output)
#'   
#' }
#' 
#' x = runModel(parInfo$guess, forcing) 
#' print(x)
#' names(x)
#' 
#' obj  = createObjectiveFunction(runModel=runModel, info=calibrationInfo, observed=observed, forcing=forcing)
#' obj2 = createObjectiveFunction(runModel=runModel, info=calibrationInfo, observed=observed, forcing=forcing, aggregate=TRUE)
#' cat("Starting calibration...\n")
#' 
#' calib = calibrate(par=parInfo$guess, fn=obj, lower=parInfo$lower, 
#'                   upper=parInfo$upper, phases=parInfo$phase,
#'                   control=list(weights=calibrationInfo$weights,
#'                                REPORT=10, trace=5))
#'} 
#' @export 
calibrarDemo = function(path=NULL, model=NULL,  ...) {
  
  if(is.null(path)) path = getwd()
  if(is.null(model)) {
    model = "default"
    warning("Using default demo 'PoissonMixedModel'")
  }
  
  output = switch(model, 
                  PoissonMixedModel = .generatePoissonMixedModel(path=path, ...),
                  PredatorPrey      = .generatePredatorPreyModel(path=path, ...),
                  .generatePoissonMixedModel(path=path, ...)  
  )
  output$value = NA
  output$time = NA
  output$counts = c('function'=NA, gradient=NA)
  class(output) = c("calibrar.demo", "calibrar.results", class(output))
  return(output)                
  
}

#' @export
#' @method print calibrar.demo
print.calibrar.demo = function(x, ...) {
  print.default(x, ...)
}

