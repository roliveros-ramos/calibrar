# wrapper for optim -------------------------------------------------------

.optim = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(par)
  # defaults for optim (taken from stats::optim)
  con = list(trace = 0, fnscale = 1, parscale = rep.int(1, npar), 
             ndeps = rep.int(0.001, npar), maxit = 100L, abstol = -Inf, 
             reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
             gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1, 
             lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  if (method == "Nelder-Mead") con$maxit = 500
  if (method == "SANN") {
    con$maxit = 10000
    con$REPORT = 100
  }
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(stats::optim(par=par, fn=fn, gr=gr, method=method, lower=lower, 
                                         upper=upper, control=control, hessian=hessian))
  
  # let's ensure there's always the hessian component listed.
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}

# wrapper for nlm ---------------------------------------------------------

.nlm = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(!is.null(control$parscale)) control$typsize     = control$parscale
  if(!is.null(control$trace))    control$print.level = min(max(0, control$trace), 2)
  if(!is.null(control$reltol))   control$gradtol     = control$reltol
  if(!is.null(control$maxit))    control$iterlim     = control$maxit
  control[c("parscale", "trace", "reltol", "maxit")] = NULL
  
  con = list(typsize = rep(1, length(par)), fscale = 1, print.level = 0, 
             ndigit = 12, gradtol = 1e-6, 
             steptol = 1e-6, iterlim = 100, check.analyticals = TRUE)
  con$stepmax = max(1000 * sqrt(sum((par/con$typsize)^2)), 1000)
  
  control = check_control(control=control, default=con)

  attr(fn, "gradient") = gr
  attr(fn, "hessian")  = NULL
  
  xoutput = suppressWarnings(stats::nlm(f=fn, p=par, hessian=hessian, 
      typsize = control$typsize, fscale=control$fscale, 
      print.level = control$print.level, ndigit = control$ndigit, 
      gradtol = control$gradtol, stepmax = control$stepmax,
      steptol = control$steptol, iterlim = control$iterlim, 
      check.analyticals = control$check.analyticals))

  xcode = xoutput$code
  msg = c("Relative gradient is close to zero, current iterate is probably solution.",
          "Successive iterates within tolerance, current iterate is probably solution.",
          "Last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.",
          "Iteration limit exceeded.",
          "maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.")
  message = msg[xcode]

  code = c(0, 0, 21, 1, 22)[xcode]
    
  output = list()
  output$par  = xoutput$estimate
  output$value = xoutput$minimum
  output$counts = c('function'=NA, gradient=NA)
  output$convergence = code
  output$message = message
  output$hessian = xoutput$hessian
  
  if(is.null(output$hessian)) output$hessian = NULL
  
  return(output)
  
}


# wrapper for nlminb ------------------------------------------------------


.nlminb = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(!is.null(control$maxit))  control$eval.max = 1.5*control$maxit
  if(!is.null(control$maxit))  control$iter.max = control$maxit
  if(!is.null(control$trace))  control$trace    = ifelse(control$trace>0, control$REPORT, 0)
  if(!is.null(control$abstol)) control$abs.tol  = ifelse(control$abstol==-Inf, 0, control$abstol)
  if(!is.null(control$reltol)) control$rel.tol  = control$reltol
  control[c("maxit", "REPORT", "abstol", "reltol")] = NULL
  
  con = list(eval.max = 200, iter.max = 150, trace = 0, abs.tol = 0, rel.tol = 1e-10,
             x.tol = 1.5e-8, xf.tol = 2.2e-14, step.min = 1, step.max = 1, 
             scale.init = NULL, diff.g = NULL, sing.tol=NULL)
  
  control = check_control(control=control, default=con)
  if(is.null(control$scale.init)) control$scale.init = NULL
  if(is.null(control$diff.g)) control$diff.g = NULL
  if(is.null(control$sing.tol)) control$sing.tol = control$rel.tol
  
  attr(fn, "gradient") = gr
  attr(fn, "hessian")  = NULL
  
  xoutput = suppressWarnings(stats::nlminb(start=par, objective=fn, gradient=gr, hessian=NULL, 
         scale = 1, control=control, lower=lower, upper=upper))
  
  output = list()
  output$par  = xoutput$par
  output$value = xoutput$objective
  output$counts = xoutput$evaluations	
  output$convergence = xoutput$convergence
  output$message = xoutput$message
  output$hessian = NULL
  
  return(output)
  
}

