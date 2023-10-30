
# wrapper for dfoptim -----------------------------------------------------

.dfoptim = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(!is.null(control$trace)) {
    control$trace = ifelse(control$trace>1, TRUE, FALSE)
    control$info = control$trace
  } 
  
  control$ncores = control$nvar = control$stochastic = control$gradient = 
    control$gr.method = NULL
   
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
  
  
  control = check_control(control=control, default=con, minimal=FALSE, verbose=FALSE)
  
  xoutput = suppressWarnings(switch(method,
                                mads = dfoptim::mads(par=par, fn=fn, lower=lower, upper=upper, scale=1, control=control),
                                hjk  = dfoptim::hjk(par=par, fn=fn, control=control),
                                hjkb = dfoptim::hjkb(par=par, fn=fn, lower=lower, upper=upper, control=control),
                                nmk  = dfoptim::nmk(par=par, fn=fn, control=control),
                                nmkb = dfoptim::nmkb(par=par, fn=fn, lower=lower, upper=upper, control=control)
  ))
  
  output = list()
  output$par  = xoutput$par
  output$value = xoutput$value
  output$counts = c('function'=xoutput$feval, 'gradient'=0)	
  output$convergence = 0
  output$message = NULL
  output$hessian = NULL
  output$trace = xoutput$iterlog
  
  return(output)
  
}

