
# Methods for calibrar.result class ---------------------------------------

print.calibrar.result = function(x, ...) {
  
  cat("Function value:", calib$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
   
}

coef.calibrar.result = function(x, ...) {
  return(x$par)
}

plot.calibrar.result = function(x, ...) {
  return(invisible(NULL))
}

summary.calibrar.result = function(x, ...) {
  x$nphases = length(x$phases)
  x$nactive = sum(x$active)
  x$npar = length(x$par)
  class(x) = "summary.calibrar.result"
  return(x)
}

print.summary.calibrar.result = function(x, ...) {
  cat(sprintf("Calibration in %d phases.\n", x$nphases))
  cat("Function value:", calib$value, "\n")
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

print.optimES.result = function(x, ...) {
  
  cat("\nFunction value:", calib$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  if(!isTRUE(x$active$flag)) cat("* Only active parameters are shown.\n")
  
}

coef.optimES.result = function(x, ...) {
  return(x$par)
}

plot.optimES.result = function(x, ...) {
  return(invisible(NULL))
}

summary.optimES.result = function(x, ...) {
  class(x) = "summary.optimES.result"
  return(x)
}

print.summary.optimES.result = function(x, ...) {
  cat("\nFunction value:", calib$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  if(!isTRUE(x$active$flag)) cat("* Only active parameters are shown.")

  cat("Partial fitness values:\n")
  print(x$partial)
  
  cat("Counts:\n")
  print(x$counts)
  
}
