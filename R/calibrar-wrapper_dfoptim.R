
# wrapper for dfoptim -----------------------------------------------------

.dfoptim = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(!is.null(control$trace)) {
    control$trace = ifelse(control$trace>1, TRUE, FALSE)
    control$info = control$trace
  } 
  
  con = switch(method,
               mads = list(tol = 1e-06, 
                           maxfeval = 10000, 
                           maximize = FALSE, 
                           pollStyle = "lite", deltaInit = 0.01, expand = 4, 
                           lineSearch = 20, seed = 1138, trace = FALSE),
               hjk  = list(tol = 1e-06, 
                           maxfeval = Inf, 
                           maximize = FALSE, 
                           target = Inf, info = FALSE),
               hjkb = list(tol = 1e-06, 
                           maxfeval = Inf, 
                           maximize = FALSE, 
                           target = Inf, info = FALSE),
               nmk  = list(tol = 1e-06, 
                           maxfeval = min(5000, max(1500, 20 * length(par)^2)), 
                           maximize = FALSE,
                           regsimp = TRUE, restarts.max = 3, trace = FALSE),
               nmkb = list(tol = 1e-06, 
                           maxfeval = min(5000, max(1500, 20 * length(par)^2)), 
                           maximize = FALSE,
                           regsimp = TRUE, restarts.max = 3, trace = FALSE)
  )
  
  
  control = check_control(control=control, default=con)
  
  out = suppressWarnings(switch(method,
                                mads = mads(par=par, fn=fn, lower=lower, upper=upper, scale=1, control=control),
                                hjkb = hjkb(par=par, fn=fn, lower=lower, upper=upper, control=control),
                                nmkb = nmkb(par=par, fn=fn, lower=lower, upper=upper, control=control)
  ))
  
  counts = c('function'=out$feval, 'gradient'=0)
  convergence = 0
  
  output = list(par=par, value=value, counts=counts, convergence=0)
  
  return(output)
  
}

