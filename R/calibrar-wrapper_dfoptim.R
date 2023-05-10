
# wrapper for dfoptim -----------------------------------------------------

.dfoptim = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  pathTmp = getwd()               # get the current path
  on.exit(setwd(pathTmp))         # back to the original path after execution
  
  parNames = names(par)
  
  control$trace = TRUE # 
  control$info = control$trace
  
  ctl = switch(method,
               mads = list(tol = 1e-06, 
                           maxfeval = 10000, 
                           maximize = FALSE, 
                           pollStyle = "lite", deltaInit = 0.01, expand = 4, 
                           lineSearch = 20, seed = 880820, trace = FALSE),
               hjkb = list(tol = 1e-06, 
                           maxfeval = Inf, 
                           maximize = FALSE, 
                           target = Inf, info = FALSE),
               nmkb = list(tol = 1e-06, 
                           maxfeval = min(5000, max(1500, 20 * length(par)^2)), 
                           maximize = FALSE,
                           regsimp = TRUE, restarts.max = 3, trace = FALSE)
  )
  
  
  control = NULL # check!
  
  out = suppressWarnings(switch(method,
                                mads = mads(par=par, fn=fn, lower=lower, upper=upper, scale=1, control=control),
                                hjkb = hjkb(par=par, fn=fn, lower=lower, upper=upper, control=control),
                                nmkb = nmkb(par=par, fn=fn, lower=lower, upper=upper, control=control)
  ))
  
  counts = c('function'=out$feval, 'gradient'=0)
  convergence = 0
  
  output = list(par=par, value=value, counts=counts, convergence=0)
  
  # tol
  # maxfeval
  # maximize
  
  # mads control  
  # trace
  # pollStyle
  # deltaInit
  # expand
  # lineSearch
  # seed
  
  # nmkb control
  # regsimp
  # restarts.max
  # trace
  
  # hjkb control
  # target
  # info
  
  return(output)
  
}

# mads(par, fn, lower=-Inf, upper=Inf, scale=1, control = list(), ...)
# hjkb(par, fn, lower = -Inf, upper = Inf, control = list(), ...)
# nmkb(par, fn, lower=-Inf, upper=Inf, control = list(), ...)