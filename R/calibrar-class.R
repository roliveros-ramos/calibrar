
# Methods for calibrar.results class ---------------------------------------

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

#' @export
coef.calibrar.results = function(object, ...) {
  return(object$par)
}

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


#' @export
plot.calibrar.results = function(x, ...) {
  return(invisible(NULL))
}

#' @title Summary for calibration results object
#' @param object Object of class calibrar.results
#' @param ... More objects of class calibrar.results as needed for comparisons.
#' @param show_par Vector of names of positions of the parameters to show in the summary. 
#' @param par.only Show only parameters in the summary, used when more than one optimization results are summarized.
#' @keywords internal
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
  object$nphases = length(object$phases)
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

# Methods for ahres.result class ----------------------------------------

#' @export
print.ahres.result = function(x, short=FALSE, ...) {
  
  cat("\nFunction value:", x$value, "\n")
  if(!isTRUE(short)) {
    cat(sprintf("Parameters (%d of %d parameters active).\n",
                length(x$active.par), length(x$par)))
    print(x=x$par, ...)    
    if(!isTRUE(x$active.flag)) cat("* Parameters not calibrated.\n")
  }
  
}

#' @export
coef.ahres.result = function(object, ...) {
  return(object$par)
}

#' @export
plot.ahres.result = function(x, ...) {
  return(invisible(NULL))
}

#' @export
summary.ahres.result = function(object, ...) {
  class(object) = "summary.ahres.result"
  return(object)
}

#' @export 
print.summary.ahres.result = function(x, ...) {
  cat("\nFunction value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  if(!isTRUE(x$active$flag)) cat("* Only active parameters are shown.")

  cat("Partial fitness values:\n")
  print(x$partial)
  
  cat("Counts:\n")
  print(x$counts)
  
}

#' @export
dim.calibrar.results = function(x) {
  length(unlist(x$par))
}
