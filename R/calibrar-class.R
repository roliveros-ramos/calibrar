# -----------------------------------------------------------------------------
# Draft roxygen documentation for result objects and S3 methods in calibrar
# Suggested target file: R/calibrar-class.R
# -----------------------------------------------------------------------------

#' Calibration result objects
#'
#' @description
#' Objects of class `calibrar.results` are returned by [calibrate()]. They store
#' the final parameter values together with the objective value, convergence
#' information, the set of calibrated parameters, and tracing information from
#' the calibration run.
#'
#' @details
#' A `calibrar.results` object represents the final output of a complete
#' calibration workflow, potentially including multiple calibration phases.
#'
#' The most useful accessors are:
#'
#' - [print()] for a compact console summary;
#' - [summary()] for a more explicit summary, including selected parameters;
#' - [coef()] to extract the fitted parameter values;
#' - [predict()] to evaluate the fitted parameter set with the underlying model
#'   function when this is available; and
#' - [dim()] to report the number of parameters in the fitted solution.
#'
#' The [predict()] method is currently only available when the objective
#' function was built with [calibration_objFn()], because it relies on metadata
#' attached to that function.
#'
#' @docType class
#' @name calibrar.results
#' @rdname calibrar.results
#' @format A list, usually with the following components:
#' \describe{
#'   \item{par}{Best parameter values found, returned with the same structure as
#'   the input `par` supplied to [calibrate()].}
#'   \item{value}{Objective-function value at `par`.}
#'   \item{counts}{Named integer vector giving the number of function and
#'   gradient evaluations recorded for the last calibration phase.}
#'   \item{convergence}{Integer convergence code returned by the optimiser used
#'   in the final phase.}
#'   \item{message}{Character string with additional information returned by the
#'   optimiser.}
#'   \item{method}{Name of the optimisation method used in the calibration.}
#'   \item{fn}{Objective function used for the calibration.}
#'   \item{active}{Logical vector indicating which parameters were active
#'   (estimated) in the full calibration problem.}
#'   \item{elapsed}{Elapsed time for the complete calibration run.}
#'   \item{trace}{A list containing tracing information from the last phase and
#'   additional run-level metadata. Depending on the control settings and the
#'   optimiser used, this may include per-phase outputs, generation traces,
#'   partial objective values, and other diagnostic information.}
#' }
#'
#' @seealso [calibrate()], [calibration_objFn()], [summary()], [coef()],
#'   [predict()]
#' @keywords classes
NULL

#' @rdname calibrar.results
#' @title Methods for `calibrar.results` objects
#' @param x A `calibrar.results` object, or a `summary.calibrar.results` object
#'   for `print.summary.calibrar.results()`.
#' @param object A `calibrar.results` object.
#' @param ... Additional arguments passed to the underlying method.
#' @export
print.calibrar.results = function(x, ...) {
  
  cat(sprintf("Optimization using '%s' algorithm.\n", x$method))
  cat("Function value:", x$value, "\n")
  cat("Status:", x$message, "\n")
  cat("Parameters:\n")
  print(x=unlist(x$par), ...)
  if(!all(x$active)) cat("* Some parameters are not calibrated.\n")
  if(!is.null(x$trace$generations)) {
    cat(sprintf("Computation: (%d generations)\n", x$trace$generations))
  } else {
    cat("Computation:\n")
  }
  print(x$counts)
  return(invisible(NULL))
}

#' @rdname calibrar.results
#' @description
#' `coef()` returns the fitted parameter values stored in a
#' `calibrar.results` object.
#' @return `coef.calibrar.results()` returns the fitted parameter values, with
#'   the same structure as the `par` component of the object.
#' @export
coef.calibrar.results = function(object, ...) {
  return(object$par)
}

#' @rdname calibrar.results
#' @param replicates Number of replicate evaluations of the fitted parameter set
#'   to perform. If `replicates = 1`, the underlying model function is evaluated
#'   once. If `replicates > 1`, repeated evaluations are combined by taking the
#'   mean across replicates for each returned output.
#' @param na.rm Logical. Passed to the internal averaging step used when
#'   `replicates > 1`.
#' @return
#' `predict.calibrar.results()` returns the model outputs produced by the fitted
#' parameter set. When `replicates > 1`, it returns the mean output across
#' replicates and stores the full replicate-level outputs in the `.replicates`
#' component of the returned object.
#' @export
predict.calibrar.results = function(object, replicates=1, 
                                    na.rm=FALSE, ...) {

  obj = object$fn
  
  if(is.null(attr(obj, which="fn"))) {
    warning("Predict is only available for functions created with 'calibration_objFn'.")
    return(invisible(NULL))
  }
  
  fn = match.fun(attr(obj, which="fn"))  
  
  if(replicates <= 1) {
    out = fn(object$par)
  } else {
    xout = list()
    for(i in seq_len(replicates)) {
      xout[[i]] = fn(object$par)
    } 
    xout$FUN = ".mycbind"
    xout = do.call(what=mapply, args=xout)
    out = lapply(xout, FUN=.myRowMean, na.rm=na.rm)
    out$.replicates = xout
  }
  return(out)
}


#' @rdname calibrar.results
#' @description
#' `plot()` is currently a placeholder method for `calibrar.results` objects.
#' No graphical summary is produced at present.
#' @return `plot.calibrar.results()` returns `NULL`, invisibly.
#' @export
plot.calibrar.results = function(x, ...) {
  return(invisible(NULL))
}

#' @rdname calibrar.results
#' @param show_par Optional character vector of parameter names or integer vector
#'   of positions indicating which parameters should be shown in the summary.
#'   By default, all fitted parameters are shown.
#' @param par.only Logical. If `TRUE`, return only the method name and selected
#'   parameter values. This is mainly useful when summarising several
#'   `calibrar.results` objects side by side.
#' @return
#' `summary.calibrar.results()` returns either:
#' 
#' - an object of class `summary.calibrar.results` when called on a single
#'   object; or
#' - a `data.frame` with one row per result object when several
#'   `calibrar.results` objects are supplied through `...`.
#'
#' The print method formats these summaries for console display.
#' @export
summary.calibrar.results = function(object, ..., show_par=NULL, par.only=FALSE) {
  oNames = as.character(match.call())[-1]
  objs = list(...)
  useDots = all(sapply(objs, FUN=inherits, what="calibrar.results"))
  if(useDots & length(objs)>0) {
    objs = c(list(object), objs)
    out = lapply(objs, FUN=.summaryCalibrarResults, pars=show_par, par.only=par.only)
    out = do.call(rbind, out)
    rownames(out) = oNames[seq_along(objs)]
    class(out) = c("summary.calibrar.results", class(out))
    return(out)
  }
  object$nphases = length(object$trace$phases)
  object$nactive = sum(object$active)
  object$npar = length(unlist(object$par))
  object$show_pars = show_par
  object$par.only = par.only
  class(object) = "summary.calibrar.results"
  return(object)
}

.summaryCalibrarResults = function(x, pars=NULL, par.only=FALSE) {
  ox = x
  xpars = unclass(unlist(x$par))
  npar = length(xpars)
  if(is.null(pars)) pars = seq_len(npar)
  if(is.character(pars)) {
    if(all(pars %in% names(x$par))) {
      xpars = xpars[pars]
      xpars = unclass(unlist(xpars))
    }
  } else {
    if(all(pars %in% seq_len(npar)))
      xpars = unlist(unclass(xpars))[pars]
  }
  
  if(isTRUE(par.only)) {
    out = data.frame(method=x$method, t(xpars))
  } else {
    out = data.frame(method=x$method, elapsed=x$elapsed, value=x$value, fn=x$counts[1], gr=x$counts[2], t(xpars))
  }
  return(out)
}

#' @rdname calibrar.results
#' @param digits Number of significant digits to use when printing summary
#'   tables or parameter values.
#' @export
print.summary.calibrar.results = function(x, digits=3, ...) {
  if(is.data.frame(x)) {
    print.data.frame(x, digits=digits, ...)
    return(invisible())
  }
  nphases = length(x$trace$phases)
  cat(sprintf("Calibration in %d %s.\n\n", nphases, 
              ifelse(nphases==1, "phase", "phases")))
  spar = unlist(unclass(x$par))
  if(!is.null(x$show_pars)) {
    cat(sprintf("Optimal parameter values (%d of %d):\n",
                length(x$show_pars),
                length(spar)))
    print(x=spar[x$show_pars], ...)
  } else {
    cat("Optimal parameter values:\n")
    print(x=spar, ...)
  }
  cat(sprintf("\n\t%d of %d parameters have been calibrated.\n\n", 
              x$nactive, x$npar))
  cat("Method         :", x$method, "\n")
  cat("Function value :", x$value, "\n")
  cat(sprintf("Counts         : fn=%d, gr=%d\n", x$counts[1], x$counts[2]))
  if(!is.null(x$partial)) {
    val =  paste(sprintf("%0.2f", x$partial), collapse=", ")
    cat("Partial fitness values: \n", val)
    print()
  }
  cat(sprintf("Message        : %s\n", x$message))
  return(invisible())
}

#' @rdname calibrar.results
#' @description
#' `dim()` returns the number of fitted parameters after unlisting the parameter
#' structure stored in the result object.
#' @return `dim.calibrar.results()` returns a single integer.
#' @export
dim.calibrar.results = function(x) {
  length(unlist(x$par))
}


#' Evolutionary-strategy result objects
#'
#' @description
#' Objects of class `ahres.results` are returned by [ahres()]. They store the
#' fitted parameters and the main diagnostics produced by the AHR-ES optimiser.
#'
#' @details
#' These objects are specific to the AHR-ES optimisation workflow. They contain
#' the final parameter vector, the aggregated objective value, convergence and
#' iteration information, and, when relevant, partial objective values and
#' tracing information.
#'
#' The main accessors are [print()], [summary()], and [coef()]. A `plot()`
#' method is currently exported as a placeholder but does not yet produce a
#' graphical summary.
#'
#' @docType class
#' @name ahres.results
#' @rdname ahres.results
#' @format A list, usually with the following components:
#' \describe{
#'   \item{par}{Best parameter values found by the AHR-ES algorithm.}
#'   \item{value}{Aggregated objective-function value at `par`.}
#'   \item{counts}{Named integer vector giving the number of function and
#'   gradient evaluations recorded by the optimiser.}
#'   \item{convergence}{Integer convergence code. `0` indicates that the
#'   stopping criterion was reached before the generation limit; `1` indicates
#'   that the maximum number of generations or function evaluations was reached.}
#'   \item{message}{Character string describing the stopping condition.}
#'   \item{generations}{Number of generations completed.}
#'   \item{partial}{Optional vector of partial objective-function values,
#'   present when the objective function returns multiple components.}
#'   \item{trace}{Optional tracing information collected during the run.}
#'   \item{active}{Information on which parameters were active in the
#'   optimisation, as stored by [ahres()].}
#' }
#'
#' @seealso [ahres()], [summary()], [coef()]
#' @keywords classes
NULL

#' @rdname ahres.results
#' @title Methods for `ahres.results` objects
#' @param x An `ahres.results` object, or a `summary.ahres.results` object for
#'   `print.summary.ahres.results()`.
#' @param object An `ahres.results` object.
#' @param ... Additional arguments passed to the underlying method.
#' @param short Logical. If `TRUE`, print a compact summary without listing the
#'   full parameter vector.
#' @export
print.ahres.results = function(x, short=FALSE, ...) {
  
  cat("\nFunction value:", x$value, "\n")
  if(!isTRUE(short)) {
    cat(sprintf("Parameters (%d out of %d parameters are currently active).\n",
                sum(x$active$isActive), length(x$par)))
    print(x=x$par, ...)    
    if(!isTRUE(x$active$activeFlag)) cat("* Parameters not calibrated.\n")
  }
  
}

#' @rdname ahres.results
#' @description
#' `coef()` returns the fitted parameter values stored in an `ahres.results`
#' object.
#' @return `coef.ahres.results()` returns the fitted parameter vector.
#' @export
coef.ahres.results = function(object, ...) {
  return(object$par)
}

#' @rdname ahres.results
#' @description
#' `plot()` is currently a placeholder method for `ahres.results` objects.
#' No graphical summary is produced at present.
#' @return `plot.ahres.results()` returns `NULL`, invisibly.
#' @export
plot.ahres.results = function(x, ...) {
  return(invisible(NULL))
}

#' @rdname ahres.results
#' @return
#' `summary.ahres.results()` returns an object of class
#' `summary.ahres.results`.
#' @export
summary.ahres.results = function(object, ...) {
  class(object) = "summary.ahres.results"
  return(object)
}

#' @rdname ahres.results
#' @export 
print.summary.ahres.results = function(x, ...) {
  cat("\nFunction value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  if(!isTRUE(x$active$flag)) cat("* Only active parameters are shown.")

  cat("Partial fitness values:\n")
  print(x$partial)
  
  cat("Counts:\n")
  print(x$counts)
  
}
