
gradient = function(fn, x, method, control, ...) {
  UseMethod("gradient")
}

gradient.default = function(fn, x, method, control=list(), ...) {
  
  if(is.null(control$parallel)) control$parallel = FALSE
  
  if(isTRUE(control$parallel)) {
    
    df = switch(method,
                central    = .grad_central_parallel(fn, x, control=control, ...),
                forward    = .grad_simple_parallel(fn, x, side=+1, control=control, ...),
                backward   = .grad_simple_parallel(fn, x, side=-1, control=control, ...),
                richardson = .grad_richardson_parallel(fn, x, control=control, ...),
                stop("Undefined method."))
    
  } else {

    df = switch(method,
                central    = .grad_central(fn, x, control=control, ...),
                forward    = .grad_simple(fn, x, side=+1, control=control, ...),
                backward   = .grad_simple(fn, x, side=-1, control=control, ...),
                richardson = .grad_richardson(fn, x, control=control, ...),
                stop("Undefined method."))
    
  }
    
  return(df)
  
}
  

# Implementation of methods for gradient computation ----------------------

.grad_simple = function(fn, x, side=1, control=list(), ...) {
  # simple method uses n+1 function evaluation, n can be parallelized,
  # so computer time is 2 functions.
  
  n = length(x)
  
  h = .get_h(x, control, method="simple")*side

  if(length(side)==1)   side = rep(side, n)
  if(length(side)!=n)   stop("Argument 'side' should have the same length as x.")
  if(any(is.na(side)))  side[is.na(side)] = 1
  if(any(abs(side)!=1)) stop("Argument 'side' should have values +1 or -1.")
  
  fx = fn(x, ...)
  
  df = rep(NA_real_, n)
  
  for(i in seq_len(n)) {
    dx = x + h*(i == seq_len(n))
    df[i] = (fn(dx, ...) - fx)/h[i]
  }
  return(df)
}


.grad_central = function(fn, x, control=list(), ...) {
  # central method uses 2*n function evaluation, n can be parallelized,
  # so computing time is 2 functions.
  
  n = length(x)

  h = .get_h(x, control, method="simple")
  
  df = rep(NA_real_, n)
  
  for(i in seq_len(n)) {
    dp = x + h*(i == seq_len(n))
    dm = x - h*(i == seq_len(n))
    df[i] = (fn(dp, ...) - fn(dm, ...))/(2*h[i])
  }
  
  return(df)
}

.grad_richardson = function(fn, x, control=list(), ...) {
  # richardson method uses 2*r*n function evaluations, n can be parallelized,
  # so computing time is 2*r. The default is r=4, so 8 function evaluations eq,
  # and it's 4 times more costly than central or simple method.
  
  n = length(x)
  
  h = .get_h(x, control, method="start")

  DF = matrix(NA_real_, nrow=r, ncol=n)
  
  for(k in seq_len(r)) {
    
    for(i in seq_len(n)) {
      
      if((k > 1) && (abs(DF[(k - 1), i]) < 1e-20)) {
        # set derivative to zero when below 1e-20, why hardcoded?
        DF[k, i] = 0
        
      } else {
        
        dp = x + h*(i == seq_len(n))
        dm = x - h*(i == seq_len(n))
        DF[k, i] = (fn(dp, ...) - fn(dm, ...))/(2*h[i])
        
        if(any(is.na(DF[k, i])))
          stop(sprintf("function returns NA at %g distance from x.", h))
        
      }
    }
    
    h = h/v
    
  }
  
  # ponderate all estimates
  for(m in seq_len(r - 1)) {
    
    ind0 = 2:(r + 1 - m)
    ind1 = 1:(r - m)
    DF = (DF[ind0, , drop = FALSE]*(4^m) - DF[ind1, , drop = FALSE])/(4^m - 1)
    
  }
  
  return(as.numeric(DF))
  
}


# Auxiliar functions ------------------------------------------------------

.get_h = function(x, control, method) {
  
  if(method=="simple") {
    args = list(eps = 1e-08, zero.tol = sqrt(.Machine$double.eps/7e-07))
    args[names(control)] = control
    eps = args$eps 
    zero.tol = args$zero.tol # threshold for an initial value to be zero (for initial h).
    isZero = 0 + (abs(x) < zero.tol)
    h = eps*abs(x)*(1-isZero) + eps*isZero
    return(h)
  }
  
  if(method=="optextras") {
    args = list(eps = 1e-08)
    args[names(control)] = control
    eps = args$eps 
    h = eps*(abs(x) + eps)
    return(h)
  }
  
  if(method=="start") {
    args = list(eps = 1e-04, d = 1e-04, r = 4, v = 2, zero.tol = sqrt(.Machine$double.eps/7e-07))
    args[names(control)] = control
    
    d = args$d # fraction of x for initial h.
    r = args$r # number of iterations, by default 4.
    v = args$v # reduction factor of h, by default 2.
    eps = args$eps # initial epsilon when x near zero.
    zero.tol = args$zero.tol # threshold for an initial value to be zero (for initial h).
    
    isZero = 0 + (abs(x) < zero.tol)
    h = d*abs(x)*(1-isZero) + eps*isZero
    return(h)
  }
  
  stop("Undefined method for 'h' calculation.")
  
}

