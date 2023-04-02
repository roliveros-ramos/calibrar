

# Rcgmin ------------------------------------------------------------------

.Rcgmin = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for Rcgmin (taken from Rcgmin::Rcgmin)
  con = list(maxit = 500, maximize = FALSE, trace = 0, eps = 1e-07, 
               dowarn = TRUE, tol = 0, checkgrad=FALSE, checkbounds=TRUE)
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(Rcgmin::Rcgmin(par=par, fn=fn, gr=gr, lower=lower, 
                                           upper=upper, control=control))
  
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
             checkgrad=FALSE, checkbounds=TRUE, keepinputpar=FALSE)
  
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

