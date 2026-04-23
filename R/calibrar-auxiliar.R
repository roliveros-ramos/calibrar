
# calibrar_demo()

#' Create small demonstration calibration projects
#'
#' Create self-contained example projects that can be used to learn the
#' `calibrar` workflow, test installation, or prototype calibration code.
#' Each demo writes a small set of files to disk and returns the information
#' needed to read the setup table, load the observed data, define the model
#' interface, and run a calibration.
#'
#' @param path Directory where the demo project should be created. If `NULL`,
#'   the current working directory is used.
#' @param model Name of the demo model to generate. Supported values currently
#'   include `"PoissonMixedModel"`, `"PredatorPrey"`, `"SIR"`, and
#'   `"IBMLotkaVolterra"`.
#' @param ... Additional model-specific arguments passed to the corresponding
#'   demo generator, for example the simulation length or the true parameter
#'   values used to generate the synthetic observations.
#'
#' @return An object of class `calibrar.demo`, which also inherits from
#'   [calibrar.results()]. The returned list contains the demo project path,
#'   the true parameter values used to generate the synthetic data, the path to
#'   the calibration setup file, suggested starting values and bounds, phase
#'   information, model-specific constants, and a precomputed objective value at
#'   the true parameter set.
#'
#' @details
#' The generated project is meant to be a compact, runnable example of the
#' `calibrar` workflow rather than a realistic scientific application. The
#' returned object typically contains enough information to proceed through the
#' standard helper functions:
#'
#' - [calibration_setup()] to read the setup table,
#' - [calibration_data()] to load the observed data,
#' - [calibration_objFn()] to build the objective function, and
#' - [calibrate()] or [optim2()] to run the optimisation.
#'
#' The exact files created depend on `model`, but they usually include a
#' calibration settings table, one or more observation files, and any required
#' forcing or master data files.
#'
#' The currently available demos are:
#'
#' \describe{
#'   \item{`"PoissonMixedModel"`}{A Poisson autoregressive mixed model for the
#'   dynamics of a population observed at several sites over time.}
#'   \item{`"PredatorPrey"`}{A deterministic Lotka--Volterra predator--prey
#'   system solved numerically.}
#'   \item{`"SIR"`}{A simple susceptible--infected--recovered epidemic model
#'   solved numerically.}
#'   \item{`"IBMLotkaVolterra"`}{A stochastic individual-based predator--prey
#'   example.}
#' }
#'
#' @seealso
#' [calibration_setup()], [calibration_data()], [calibration_objFn()],
#' [calibrate()], [summary.calibrar.results()], [print.calibrar.results()]
#'
#' @examples
#' tmp = tempdir()
#' demo = calibrar_demo(path = tmp, model = "PoissonMixedModel", T = 8, L = 3)
#' print(demo)
#'
#' setup = calibration_setup(demo$setup)
#' observed = calibration_data(setup, path = demo$path)
#'
#' @export
calibrar_demo = function(path=NULL, model=NULL, ...) {
  
  if(is.null(path)) path = getwd()
  if(is.null(model)) {
    model = "PoissonMixedModel"
    warning("Using default demo 'PoissonMixedModel'")
  }
  
  output = switch(model, 
                  PoissonMixedModel = .generatePoissonMixedModel(path=path, ...),
                  PredatorPrey      = .generatePredatorPreyModel(path=path, ...),
                  SIR               = .generateSIRModel(path=path, ...),
                  IBMLotkaVolterra  = .generateIBMLotkaVolterra(path, ...),
                  stop(sprintf("Model '%s' is not defined.", model))
                 )
 
  output$elapsed = NA
  output$counts = c('function'=NA, gradient=NA)
  output$method = "data"
  output$model  = model
  class(output) = c("calibrar.demo", "calibrar.results", class(output))
  return(output)                
  
}

#' @rdname calibrar_demo
#' @param x A `calibrar.demo` object returned by [calibrar_demo()].
#' @param ... Further arguments passed to [print()].
#' @export
#' @method print calibrar.demo
print.calibrar.demo = function(x, ...) {

  cat(sprintf("Calibration demo using model '%s'.\n", x$model))
  if(!is.na(x$value)) cat("Target value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=unlist(x$par), ...)
  cat("Lower bounds:\n")
  print(x=unlist(x$lower), ...)
  cat("Upper bounds:\n")
  print(x=unlist(x$upper), ...)
  cat("Phases:\n")
  print(x=unlist(x$phase), ...)

  cat("Guess (starting point):\n")
  print(x=unlist(x$guess), ...)
  
  cat(sprintf("\nFiles location: %s\n", x$path))
  
  return(invisible())
  
}



# Some other useful stuff -------------------------------------------------

# .get_command_argument()

#' Extract a named command-line argument
#'
#' Parse command-line arguments and return the value associated with a named
#' option. This helper is primarily intended for scripting workflows that call
#' R from the shell and pass arguments such as `--config=file.txt`.
#'
#' @param x Character vector of command-line arguments, typically obtained from
#'   `commandArgs()`.
#' @param argument Name of the argument to retrieve, without the prefix.
#' @param prefix Prefix used for arguments of interest. The default is `"--"`.
#' @param default Value returned when the argument is not found. The default is
#'   `FALSE`.
#' @param verbose Logical. If `TRUE`, a warning is issued when the requested
#'   argument is not found.
#'
#' @return
#' The parsed argument value. If the argument is present without an explicit
#' `=value` part, `TRUE` is returned. If the argument is absent, `default` is
#' returned. Parsed values are passed through an internal type guesser, so
#' logical and numeric values may be returned with their corresponding R types.
#'
#' @details
#' The function looks for entries matching either `prefix + argument` or
#' `prefix + argument = value`. For example, with `prefix = "--"` and
#' `argument = "config"`, both `--config` and `--config=file.txt` are
#' recognised.
#'
#' This function is exported because it can be useful in external packages or
#' scripts that build command-line interfaces around `calibrar` workflows.
#'
#' @seealso
#' [.read_configuration()]
#'
#' @examples
#' .get_command_argument(commandArgs(), "interactive")
#' .get_command_argument(commandArgs(), "RStudio")
#' .get_command_argument(commandArgs(), "RStudio", prefix = "")
#' .get_command_argument(commandArgs(), "vanilla")
#' .get_command_argument("--control.file=baz.txt", "control.file")
#'
#' @export
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

# .read_configuration()

#' Read a simple key-value configuration file
#'
#' Read a plain-text configuration file in which each non-empty, non-comment
#' line defines a key and one or more values separated by a recognised
#' delimiter. Parsed values are returned as a named list.
#'
#' @param file Path to the configuration file.
#' @param recursive Logical. If `TRUE`, keys matching `conf.key` are treated as
#'   references to additional configuration files and are read recursively.
#' @param keep.names Logical. If `TRUE`, keys are kept as written. If `FALSE`,
#'   they are converted to lower case.
#' @param conf.key Optional pattern identifying keys that should be interpreted
#'   as additional configuration files.
#' @param ... Reserved for future extensions.
#'
#' @return
#' A named list containing the parsed configuration values. Scalars are returned
#' as length-one vectors where appropriate, and simple values are converted with
#' `type.convert()`. If `file` does not exist, `NULL` is returned with a
#' warning.
#'
#' @details
#' The reader accepts simple lines of the form `key separator value`, where the
#' separator may vary across lines. The current parser guesses the separator
#' among common choices such as `=`, `;`, `,`, `:`, and tab.
#'
#' Lines beginning with punctuation after trimming are ignored, and trailing
#' comments are removed. When `conf.key` is supplied, matching entries are
#' treated as file references. Their paths are resolved relative to the parent
#' configuration file, and the referenced files can be merged recursively.
#'
#' This helper is intentionally lightweight. It is useful for small calibration
#' or script configuration files, but it is not intended to replace richer
#' configuration formats such as YAML.
#'
#' @seealso
#' [.get_command_argument()]
#'
#' @examples
#' tf = tempfile(fileext = ".txt")
#' writeLines(c(
#'   "ncores = 4",
#'   "verbose = TRUE",
#'   "method = AHR-ES"
#' ), tf)
#'
#' cfg = .read_configuration(tf)
#' str(cfg)
#'
#' @export
.read_configuration = function(file, recursive=TRUE, keep.names = TRUE, conf.key=NULL, ...) {
  
  if(!is.null(attr(file, "path"))) file = c(file.path(attr(file, "path"), file))
  if(!file.exists(file)) {
    warning(sprintf("configuration file '%s' does not exist.", file))
    return(NULL)
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
    ckey = unique(names(values[xpos]))
    msg = sprintf("Duplicated key value for '%s', using first ocurrence (%s).", ckey, xpos[1])
    if(length(xpos)>1) {
      xpos = xpos[1]
      warning(msg)
    }
      
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
  
  # check for duplicates
  ind = duplicated(names(values))
  if(sum(ind) > 0) {
    values = values[!ind]
    dup = paste(names(config)[ind], collapse="\n")
    msg = sprintf("Removing %d duplicated values:\n %s", sum(ind), dup)
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


.guessSeparator = function(Line, nocolon=FALSE) {
  SEPARATORS = c(equal = "=", semicolon = ";",
                 coma = ",", colon = ":", tab = "\t")
  if(isTRUE(nocolon)) SEPARATORS = SEPARATORS[-4]
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
  valueseparator = .guessSeparator(values, nocolon=TRUE)
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
  if(length(x) != 1) return(x)
  xind0 = file.exists(file.path(path, x)) 
  xind1 = any(file.exists(dir(path=file.path(path, dirname(x)), 
                              pattern=sprintf("^%s\\.", basename(x)), full.names=TRUE)))
  if(xind0 || xind1 || isTRUE(force)) 
    attr(x, "path") = normalizePath(path, winslash = "/", mustWork = FALSE)
  return(x)
}
