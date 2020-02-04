
gradient = function(fn, x, method, control, ...) {
  UseMethod("gradient")
}

gradient.default = function(fn, x, method, control, ...) {
  
  df = switch(method,
              central    = .grad_central(fn, x, control=list(), ...),
              forward    = .grad_simple(fn, x, side=+1, control=list(), ...),
              backward   = .grad_simple(fn, x, side=-1, control=list(), ...),
              richardson = .grad_richardson(fn, x, control=list(), ...),
              stop("Undefined method."))
  
  return(df)
  
}
  


# Implementation of methods for gradient computation ----------------------

.grad_simple = function(fn, x, side=1, control=list(), ...) {
  
  n = length(x)
  
  args = list(eps = 1e-08)
  args[names(control)] = control
  eps = args$eps # initial epsilon when x near zero.

  if(length(side)==1)   side = rep(side, n)
  if(length(side)!=n)   stop("Argument 'side' should have the same length as x.")
  if(any(is.na(side)))  side[is.na(side)] = 1
  if(any(abs(side)!=1)) stop("Argument 'side' should have values +1 or -1.")
  
  fx = fn(x, ...)
  
  df = rep(NA_real_, n)
  h = (eps * (abs(x) + eps))*side
  
  for(i in seq_len(n)) {
    dx = x + h*(i == seq_len(n))
    df[i] = (fn(dx, ...) - fx)/h[i]
  }
  return(df)
}


.grad_central = function(fn, x, control=list(), ...) {
  
  n = length(x)

  args = list(eps = 1e-08)
  args[names(control)] = control
  eps = args$eps # initial epsilon when x near zero.
  
  df = rep(NA_real_, n)
  h = 0.5 * eps * (abs(x) + 1)
  
  for(i in seq_len(n)) {
    dp = x + h*(i == seq_len(n))
    dm = x - h*(i == seq_len(n))
    df[i] = (fn(dp, ...) - fn(dm, ...))/(2*h[i])
  }
  
  return(df)
}

.grad_richardson = function(fn, x, control=list(), ...) {
  
  n = length(x)
  
  args = list(eps = 1e-04, d = 1e-04, r = 4, v = 2, zero.tol = sqrt(.Machine$double.eps/7e-07))
  args[names(control)] = control
  
  d = args$d # fraction of x for initial h.
  r = args$r # number of iterations, by default 4.
  v = args$v # reduction factor of h, by default 2.
  eps = args$eps # initial epsilon when x near zero.
  zero.tol = args$zero.tol # threshold for an initial value to be zero (for initial h).

  isZero = 0 + (abs(x) < zero.tol)
  h = d*abs(x)*(1-isZero) + eps*isZero

  A = matrix(NA_real_, nrow=r, ncol=n)
  
  for(k in seq_len(r)) {
    
    for(i in seq_len(n)) {
      
      if((k > 1) && (abs(A[(k - 1), i]) < 1e-20)) {
        # set derivative to zero when below 1e-20, why hardcoded?
        A[k, i] = 0
        
      } else {
        
        dp = x + h*(i == seq_len(n))
        dm = x - h*(i == seq_len(n))
        A[k, i] = (fn(dp, ...) - fn(dm, ...))/(2*h[i])
        
        if(any(is.na(A[k, i])))
          stop(sprintf("function returns NA at %g distance from x.", h))
        
      }
    }
    
    h = h/v
    
  }
  
  # ponderate all estimates
  for(m in seq_len(r - 1)) {
    
    ind0 = 2:(r + 1 - m)
    ind1 = 1:(r - m)
    A = (A[ind0, , drop = FALSE]*(4^m) - A[ind1, , drop = FALSE])/(4^m - 1)
    
  }
  
  return(as.numeric(A))
  
}


