
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
#' set.seed(880820)
#' path = NULL # NULL to use the current directory
#' # create the demonstration files
#' demo = calibrar_demo(model="PoissonMixedModel", L=5, T=100) 
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
#' @export 
calibrar_demo = function(path=NULL, model=NULL,  ...) {
  
  if(is.null(path)) path = getwd()
  if(is.null(model)) {
    model = "default"
    warning("Using default demo 'PoissonMixedModel'")
  }
  
  output = switch(model, 
                  PoissonMixedModel = .generatePoissonMixedModel(path=path, ...),
                  PredatorPrey      = .generatePredatorPreyModel(path=path, ...),
                  IBMLotkaVolterra  = .generateIBMLotkaVolterra(path, ...),
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



# Some other useful stuff -------------------------------------------------

#' Get an specific argument from the command line
#'
#' @param x The command line arguments, from \code{x = commandArgs()}
#' @param argument The name of the argument.
#' @param prefix The prefix to any argument of interest, the default is "--"
#' @param default Default value to return is argument is missing, default to FALSE.
#' @param verbose Boolean, if TRUE, shows a warning when the parameter is not found.
#' 
#' @return  The value of the argument, assumed to be followed after '=' or, TRUE if nothing but the argument was found. 
#' If the argument is not found, FALSE is returned.
#' @export
#'
#' @examples 
#' .get_command_argument(commandArgs(), "interactive")
#' .get_command_argument(commandArgs(), "RStudio")
#' .get_command_argument(commandArgs(), "RStudio", prefix="")
#' .get_command_argument(commandArgs(), "vanilla")
#' .get_command_argument("--control.file=baz.txt", "control.file")
.get_command_argument = function(x, argument, prefix="--", default=FALSE, verbose=FALSE) {
  
  xdefault = if(is.null(default)) "NULL" else as.character(default)
  
  w_nomatch = sprintf("No argument matches '%s', returning %s.", argument, xdefault)
  e_dumatch = sprintf("More than one arguments match '%s', provide the exact name or check for duplicates.", argument)
  
  pctl0 = sprintf("%s%s", prefix, argument)
  pctl1 = sprintf("%s%s=", prefix, argument)
  ind0 = grepl(x=x, pattern=pctl0)
  ind1 = grepl(x=x, pattern=pctl1)
  if(sum(ind0)==0) {
    if(isTRUE(verbose)) warning(w_nomatch)
    return(default)
  }
  if(sum(ind1)==0) {
    # exist but not with assigned value
    ctl = x[ind0]
    ctl = gsub(x=ctl, pattern=pctl0, replacement="")
    ind = which(nchar(ctl)==0)
    if(length(ind)==0) {
      if(isTRUE(verbose)) warning(w_nomatch)
      return(default)
    }
    if(length(ind)!=1) stop(e_dumatch)
    return(TRUE)
  } 
  if(sum(ind1)!=1) stop(e_dumatch)
  ctl = x[which(ind1)]
  ctl = gsub(x=ctl, pattern=pctl1, replacement="")
  ctl = .guessType(ctl)
  return(ctl)
}

#' Read a configuration file.
#' 
#' File is expected to have lines of the form 'key SEP value' where key is the
#' name of the parameter, SEP a separator (can be '=' ',', ';') and value the value
#' of the parameter itself. The SEP for each line is determined and parameters values are
#' returned as a list. 
#' 
#' @param file File to read the configuration
#'
#' @param recursive Should 'osmose.configuration' keys be read as additional configuration files? Default is TRUE.
#' @param keep.names Should names be kept as they are? By default, are converted to lower case, as is expected in OSMOSE.
#' @param ... Additional arguments, not currently in use.
#'
#' @export
.read_configuration = function(file, recursive=TRUE, keep.names = TRUE, conf.key=NULL, ...) {
  
  if(!is.null(attr(file, "path"))) file = c(file.path(attr(file, "path"), file))
  if(!file.exists(file)) {
    warning(sprintf("configuration file '%s' does not exist.", file))
    return(NULL)
  }
  
  .guessSeparator = function(Line){
    SEPARATORS = c(equal = "=", semicolon = ";",
                   coma = ",", colon = ":", tab = "\t")
    guess = which.min(nchar(lapply(str_split(Line,SEPARATORS), "[", i = 1)))
    separator = SEPARATORS[guess]
    
    return(separator)
  }
  
  .getKey = function(Line, KeySeparator) {
    Key = str_split(Line, KeySeparator)[[1]][1]
    return(str_trim(Key))
  }
  
  .getValues = function(x, KeySeparator){
    start = str_locate(x, pattern=KeySeparator)[1,1]
    if(is.na(start)) return(NULL)
    values = stringr::str_sub(x, start+1, nchar(x))
    valueseparator = .guessSeparator(values)
    values = stringr::str_trim(str_split(values, valueseparator)[[1]])
    values = values[nchar(values)!=0]
    values = .guessType(values)
    return(values)
  }
  
  .comment_trim = function(x, char="#") {
    start = str_locate(x, pattern=char)[1,1]
    if(is.na(start)) return(x)
    return(str_sub(x, 1, start - 1))
  }
  
  .addPath = function(x, path, force=FALSE) {
    if(is.null(x)) return(x)
    if(!is.character(x)) return(x)
    # if(!is.null(attr(x, "path"))) 
    # path = file.path(attr(x, "path"), path)
    if(file.exists(file.path(path, x)) | isTRUE(force)) 
      attr(x, "path") = normalizePath(path, winslash = "/", mustWork = FALSE)
    return(x)
  }
  
  
  config = readLines(file) # read lines
  config = lapply(config, .comment_trim) # remove comments
  config = lapply(config, str_trim)
  config[grep("^[[:punct:]]", config)] = NULL
  config = config[nchar(config)!=0]
  
  keySeparator  = sapply(config, .guessSeparator)
  key           = mapply(.getKey, config, keySeparator)
  values        = mapply(.getValues, config, keySeparator, SIMPLIFY = FALSE)
  
  names(values) = if(isTRUE(keep.names)) key else tolower(key)
  
  if(is.null(conf.key)) return(values)
  
  force = grepl(names(values), pattern=conf.key)
  values = mapply(FUN=.addPath, x=values, force=force, 
                  MoreArgs=list(path = dirname(file)), SIMPLIFY = FALSE)
  
  ii = grep(names(values), pattern=conf.key)
  if(length(ii) == 0 | !isTRUE(recursive)) return(values)
  
  ovalues = values
  
  while(length(ii) > 0) {
    
    iifile = ovalues[[ii[1]]]
    xpos = grep(names(values), pattern=names(ovalues)[ii[1]])
    ifile = file.path(attr(iifile, "path"), iifile)
    if(!file.exists(ifile)) {
      msg = sprintf("Configuration file '%s' not found. \n File '%s' not found in %s.", 
                    names(ovalues)[ii[1]], iifile, attr(iifile, "path"))
      warning(msg, immediate. = TRUE)
      message(sprintf("Skipping read of file '%s' (not found).", iifile))
      xval = NULL
    } else {
      xval = .read_configuration(ifile, recursive=TRUE, keep.names = keep.names, conf.key=conf.key, ...)
    }
    values = append(values, xval, xpos)
    ii = ii[-1]
  }
  
  return(values)
}

.guessType = function(x) {
  xx = tolower(x)
  if(identical(xx, "null")) return(NULL)
  if(identical(xx, "na")) return(NULL)
  x[xx=="true"] = "TRUE"
  x[xx=="false"] = "FALSE"
  x[xx=="na"] = "NA"
  x = x[xx!="null"]
  out = type.convert(c(x), as.is = TRUE)
  return(out)
}

