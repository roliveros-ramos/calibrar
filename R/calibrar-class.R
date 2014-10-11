
# Methods for calibrar.result class ---------------------------------------

print.calibrar.result = function(x, ...) {
  
  cat("Function value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  cat("* not calibrated parameters")
   
}

coef.calibrar.result = function(object, ...) {
  return(object$par)
}

plot.calibrar.result = function(x, ...) {
  return(invisible(NULL))
}

summary.calibrar.result = function(object, ...) {
  object$nphases = length(object$phases)
  object$nactive = sum(object$active)
  object$npar = length(object$par)
  class(object) = "summary.calibrar.result"
  return(object)
}

print.summary.calibrar.result = function(x, ...) {
  cat(sprintf("Calibration in %d phases.\n", x$nphases))
  cat("Function value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  cat(sprintf("\n\t%d of %d parameters have been calibrated.\n\n", 
              x$nactive, x$npar))
  cat("Counts:\n")
  print(x$counts)
  cat("Partial fitness values:\n")
  print(x$partial)
}

# Methods for optimES.result class ----------------------------------------

print.optimES.result = function(x, short=FALSE, ...) {
  
  cat("\nFunction value:", x$value, "\n")
  if(!isTRUE(short)) {
    cat(sprintf("Parameters (%d of %d parameters active).\n",
                length(x$active$par), length(x$par)))
    print(x=x$par, ...)    
    if(!isTRUE(x$active$flag)) cat("* Parameters not calibrated.\n")
  }
  
}

coef.optimES.result = function(object, ...) {
  return(object$par)
}

plot.optimES.result = function(x, ...) {
  return(invisible(NULL))
}

summary.optimES.result = function(object, ...) {
  class(object) = "summary.optimES.result"
  return(object)
}

print.summary.optimES.result = function(x, ...) {
  cat("\nFunction value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  if(!isTRUE(x$active$flag)) cat("* Only active parameters are shown.")

  cat("Partial fitness values:\n")
  print(x$partial)
  
  cat("Counts:\n")
  print(x$counts)
  
}
