

# Rcgmin ------------------------------------------------------------------

.Rcgmin = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for Rcgmin (taken from Rcgmin::Rcgmin)
  con = list(maxit = 500, maximize = FALSE, trace = 0, eps = 1e-07, 
               dowarn = TRUE, tol = 0, checkgrad=FALSE)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(optimx::Rcgmin(par=par, fn=fn, gr=gr, lower=lower, 
                                           upper=upper, control=control))
  names(output$counts) = c('function', 'gradient')
  # let's ensure there's always the hessian component listed.
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}

.Rcgmin_old = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for Rcgmin (taken from Rcgmin::Rcgmin)
  con = list(maxit = 500, maximize = FALSE, trace = 0, eps = 1e-07, 
             dowarn = TRUE, tol = 0, checkgrad=FALSE)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(Rcgmin::Rcgmin(par=par, fn=fn, gr=gr, lower=lower, 
                                           upper=upper, control=control))
  names(output$counts) = c('function', 'gradient')
  # let's ensure there's always the hessian component listed.
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}


# Rvmmin ------------------------------------------------------------------

.Rvmmin = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for Rvmmin (taken from Rvmmin::Rvmmin)
  con = list(maxit = 500 + 2L*npar, maxfeval = 3000 + 10L*npar, maximize = FALSE, 
             trace = 0, eps = 1e-07, dowarn = TRUE, acctol = 1e-04, 
             stepredn = 0.2, reltest = 100, stopbadupdate = TRUE, 
             checkgrad=FALSE)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(Rvmmin(par=par, fn=fn, gr=gr, lower=lower, 
                                           upper=upper, control=control))
  
  # let's ensure there's always the hessian component listed.
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}

.Rvmmin_ror = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for Rvmmin (taken from Rvmmin::Rvmmin)
  con = list(maxit = 500 + 2L*npar, maxfeval = 3000 + 10L*npar, maximize = FALSE, 
             trace = 0, eps = 1e-07, dowarn = TRUE, acctol = 1e-04, 
             stepredn = 0.2, reltest = 100, stopbadupdate = TRUE, 
             checkgrad=FALSE)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(Rvmmin_ror(par=par, fn=fn, gr=gr, lower=lower, 
                                    upper=upper, control=control))
  
  # let's ensure there's always the hessian component listed.
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}

.Rvmmin_old = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for Rvmmin (taken from Rvmmin::Rvmmin)
  con = list(maxit = 500 + 2L*npar, maxfeval = 3000 + 10L*npar, maximize = FALSE, 
             trace = 0, eps = 1e-07, dowarn = TRUE, acctol = 1e-04, 
             stepredn = 0.2, reltest = 100, stopbadupdate = TRUE, 
             checkgrad=FALSE)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(Rvmmin::Rvmmin(par=par, fn=fn, gr=gr, lower=lower, 
                                       upper=upper, control=control))
  
  # let's ensure there's always the hessian component listed.
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}

# optimr::hjn -------------------------------------------------------------

.hjn = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for Rvmmin (taken from optimr::hjn)
  con = list(trace=0, stepsize=1, stepredn=0.1, maxfeval=2000*npar,
             eps = 1e-07)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(optimr::hjn(par=par, fn=fn, lower=lower, 
                                        upper=upper, control=control))
  
  output$message = NA
  # let's ensure there's always the hessian component listed.
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}


# BB:spg ------------------------------------------------------------------

.spg = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  
  quiet = TRUE
  if(!is.null(control$trace))  quiet = ifelse(control$trace>0, TRUE, FALSE)
  if(!is.null(control$trace))  control$trace = ifelse(control$trace>1, TRUE, FALSE)
  
  method = if(!is.null(control$spg.method)) control$spg.method else 3
  
  # defaults for spg (taken from BB::spg)
  con = list(M = 10, maxit = 1500, ftol = 1e-10, gtol = 1e-05, 
             maxfeval = 10000, maximize = FALSE, trace = FALSE, triter = 10, 
             eps = 1e-07, checkGrad = NULL, checkGrad.tol = 1e-06)
  
  control = check_control(control=control, default=con)
  
  xoutput = suppressWarnings(BB::spg(par=par, fn=fn, gr=gr, lower=lower, 
                                    upper=upper, control=control, 
                                    method=method, quiet=quiet))
  
  output = list()
  output$par  = xoutput$par
  output$value = xoutput$value
  output$counts = c('function'=xoutput$feval, 'gradient'=NA)	
  output$convergence = xoutput$convergence
  output$message = xoutput$message
  output$hessian = NULL
  
  return(output)
  
}


# wrapper for lbfgsb3 -----------------------------------------------------

.lbfgsb3 = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  con = list(trace = 0L, maxit = 100L, iprint = -1L, lmm = 5, 
             factr = 1e+07, pgtol = 0, reltol = 0, abstol = 0, info = FALSE)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(lbfgsb3c::lbfgsb3(par=par, fn=fn, gr=gr, lower=lower,
                                               upper=upper, control=control))
  
  return(output)
  
}
