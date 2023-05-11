#' Numerical computation of the gradient
#' 
#' This function calculates the gradient of a function, numerically, including the possibility
#' of doing it in parallel.
#'
#' @param fn The function
#' @param x The value to compute the gradient at.
#' @param method The method used. Currently implemented: central, backward, forward and Richardson. See details.
#' @param control A list of control arguments. Notably, control$parallel=TRUE activates parallel computation.
#' @param parallel Boolean, should numerical derivatives be calculated in parallel?
#' @param ... Additional arguments to be passed to \code{fn}.
#'
#' @return The gradient of \code{fn} at \code{x}.
#' @export
#'
#' @examples gradient(fn=function(x) sum(x^3))
gradient = function(fn, x, method, control, parallel, ...) {
  UseMethod("gradient")
}

#' @export
gradient.default = function(fn, x, method=NULL, control=list(), parallel=FALSE, ...) {

  # cluster structure must be defined outside
  # by default, not running in parallel  

  if(is.null(method)) method = "richardson"
    
    df = switch(method,
                central    = .grad_central(fn, x, control=control, parallel=parallel, ...),
                forward    = .grad_simple(fn, x, side=+1, control=control, parallel=parallel, ...),
                backward   = .grad_simple(fn, x, side=-1, control=control, parallel=parallel, ...),
                richardson = .grad_richardson(fn, x, control=control, parallel=parallel, ...),
                stop("Undefined method for gradient computation."))
    

  return(df)
  
}

# Implementation of methods for gradient computation ----------------------

.grad_simple = function(fn, x, side=1, control=list(), parallel=FALSE, ...) {
  # simple method uses n+1 function evaluation, n+1 can be parallelized,
  # so computer time is 1 function evaluation.
  
  n = length(x)
  
  if(length(side)==1)   side = rep(side, n)
  if(length(side)!=n)   stop("Argument 'side' should have the same length as x.")
  if(any(is.na(side)))  side[is.na(side)] = 1
  if(any(abs(side)!=1)) stop("Argument 'side' should have values +1 or -1.")
  
  h = .get_h(x, control, method="simple")*side
  
  if(!isTRUE(parallel)) {
    
    fx = fn(x, ...)
    df = rep(NA_real_, n)
    for(i in seq_len(n)) {
      dx = h*(i == seq_len(n))
      df[i] = (fn(x + dx, ...) - fx)/h[i]
    }
    
  } else {
    
    DX = diag(n+1)
    diag(DX) = c(0, h)
    
    f = foreach(i=seq(from=1, to=n+1), .combine=c, .verbose=FALSE, .inorder=TRUE) %dopar% {
      dx  = DX[i, ]
      fi = fn(x + dx, ...)
      fi
    }
    fx = f[1] # fx = f(i=0)
    df = (f[-1] - fx)/h
    
  }

  
  return(df)
}


.grad_central = function(fn, x, control=list(), parallel=FALSE, ...) {
  # central method uses 2*n function evaluation, 2*n can be parallelized,
  # so computing time is 1 function evaluation.
  
  n = length(x)
  h = .get_h(x, control, method="simple")
  
  if(!isTRUE(parallel)) {
    
    df = rep(NA_real_, n)
    for(i in seq_len(n)) {
      dx = h*(i == seq_len(n))
      df[i] = (fn(x + dx, ...) - fn(x - dx, ...))/(2*h[i])
    }    
    
  } else {
    
    DX = diag(2*n)
    diag(DX) = c(h, -h)

    f = foreach(i=seq(from=1, to=2*n), .combine=c, .verbose=FALSE, .inorder=TRUE) %dopar% {
      dx = DX[i, ]
      fi = fn(x + dx, ...)
      fi
    }
    fp = head(f, n)
    fm = tail(f, n)
    df = (fp - fm)/(2*h)
    
  }
  
  return(df)
}


.grad_richardson = function(fn, x, control=list(), parallel=FALSE, ...) {
  # richardson method uses 2*r*n function evaluations, 2*r*n can be parallelized,
  # so computing time is 1 function evaluation.

  r = if(is.null(control$r)) 4 else control$r # number of iterations, by default 4.
  v = if(is.null(control$v)) 2 else control$v # reduction factor of h, by default 2.
  
  n = length(x)
  h = .get_h(x, control, method="start")
  
  if(!isTRUE(parallel)) {
    
    DF = matrix(NA_real_, nrow=r, ncol=n)
    for(k in seq_len(r)) {
      # gradient as usual
      for(i in seq_len(n)) {
        if((k > 1) && (abs(DF[(k - 1), i]) < 1e-20)) {
          # set derivative to zero when below 1e-20, why hardcoded?
          DF[k, i] = 0
        } else {
          dx = h*(i == seq_len(n))
          DF[k, i] = (fn(x + dx, ...) - fn(x - dx, ...))/(2*h[i])
          if(any(is.na(DF[k, i])))
            stop(sprintf("function returns NA at %g distance from x.", h))
        }
      }
      # reduce h and try again
      h = h/v
    }
    
  } else {
    
    H = matrix(NA_real_, nrow=r, ncol=n)
    for(k in seq_len(r)) H[k, ] = h/(v^(k-1))
    
    DX = NULL
    for(k in seq_len(r)) {
      HX = diag(n)
      diag(HX) = H[k, ]
      DX = rbind(DX, HX, -HX)
    }
      
    f = foreach(i=seq(from=1, to=2*r*n), .combine=c, .verbose=FALSE, .inorder=TRUE) %dopar% {
      dx  = DX[i, ]
      fi = fn(x + dx, ...)
      fi
    }
    
    f = matrix(f, nrow=r, ncol=2*n, byrow = TRUE)
    fp = f[, 1:n, drop=FALSE]
    fm = f[, -(1:n), drop=FALSE]
    DF = (fp - fm)/(2*H)
    
  } # end of parallel conditional.
  
  # ponderate all estimates (same for all cases)
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


# Notes for future self ---------------------------------------------------

# args: skeleton, active 
# # here we have to transform fn so it takes .i and change folder to run if necessary.
# skeleton = skeleton
# if(is.null(skeleton)) skeleton = as.relistable(par)
# 
# npar = length(par)
# 
# # check active parameters
# active = .checkActive(active=active, npar=npar)
# isActive = which(active)
# activeFlag = isTRUE(all(active))
# 
# # update to active parameters only
# guess  = par
# par    = guess[isActive]
# 
# npar = length(par)
# 
# force(replicates)
# 
# # closure for function evaluation
# fn   = match.fun(fn)
# 
# # here we modify f so:
# # 1. use 'isActive' to mask some parameters ('guess' is reference)
# # 2. is re-listed according to skeleton
# # 3. is re-scaled according to control$fnscale
# # 4. is evaluated 'replicates' times
# # 5. is evaluated in the 'control$run/.i' folder.
# if(is.null(control$master)) {
#   
#   fn1  = function(x, .i=0) {
#     
#     parx = guess
#     parx[isActive] = par
#     parx = relist(flesh = parx, skeleton = skeleton)
#     output = NULL
#     for(i in seq_len(replicates)) {
#       out = fn(parx, ...)/control$fnscale
#       output = rbind(output, t(as.matrix(out)))
#     }
#     return(as.numeric(colMeans(output)))
#     
#   }
#   
# } else {
#   
#   fn1  = function(par, .i=0) {
#     
#     cwd = getwd()
#     on.exit(setwd(cwd))
#     .setWorkDir(control$run, i=.i) 
#     
#     parx = guess
#     parx[isActive] = par
#     parx = relist(flesh = parx, skeleton = skeleton)
#     output = NULL
#     for(i in seq_len(replicates)) {
#       out = fn(parx, ...)/control$fnscale
#       output = rbind(output, t(as.matrix(out)))
#     }
#     return(as.numeric(colMeans(output)))
#     
#   }
#   
# } # end fn1
# 
# # here master folder is used for initialization. Maybe it's already initialized (from higher level call),
# # but it may not.
# # we return the gradient in the real dimension, setting to zero non-active components
