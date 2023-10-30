
.continueEvolution = function(opt, control) {
  
  if(opt$gen == 0) return(TRUE) 
  
  termination = control$termination
  if(is.null(termination)) termination = 0
  tt1 = Sys.time()
  out = switch(as.character(termination), 
               "0" = TRUE,
               "1" = (opt$step >= control$convergence),
               "2" = !.N_stop(opt$sstop[1:opt$gen], N=control$max_no_improvement),
               "3" = !.N_stop(opt$sstop[1:opt$gen], N=control$max_no_improvement),
               "4" = !.N_stop(opt$sstop[1:opt$gen], N=control$fn_smoothing),
               stop("Invalid termination criteria selected.")
  )
  out = (opt$gen <= (control$maxgen - 1)) & out
  tt2 = Sys.time()
  attr(out, "elapsed") = format_difftime(tt1, tt2, value=TRUE)
  return(out)
}

# Internal stop criteria --------------------------------------------------

smooth_stop = function(opt, control) {
  
  x      = opt$the_values[1:opt$gen]
  reltol = control$reltol
  
  out = switch(as.character(control$termination), 
               "0" = FALSE,
               "1" = FALSE,
               "2" = .smooth_stop2(x=x, reltol=reltol, N=control$fn_smoothing),
               "3" = .smooth_stop3(x=x, reltol=reltol, N=control$fn_smoothing),
               "4" = .smooth_stop4(x=x, reltol=reltol, N=control$fn_smoothing),
               stop("Invalid termination criteria selected.")
  )
  return(out)
}

.the_stop = function(gen, x, method, N) {
  
  x      = x[1:gen]
  reltol = sqrt(.Machine$double.eps)
  
  out = switch(as.character(method), 
               "0" = FALSE,
               "1" = FALSE,
               "2" = .smooth_stop2(x=x, reltol=reltol, N=N),
               "3" = .smooth_stop3(x=x, reltol=reltol, N=N),
               "4" = .smooth_stop4(x=x, reltol=reltol, N=N),
               stop("Invalid termination criteria selected.")
  )
  return(out)
}

the_stop = function(x, method, N) {
  sapply(seq_along(x), FUN=.the_stop, x=x, method=method, N=N)
}

.smooth_stop2 = function(x, reltol=sqrt(.Machine$double.eps), N=10) {
  
  if(length(x) < 3*N) return(FALSE)
  
  x = .tail(x, 3*N)
  
  ind = which(is.infinite(x))
  if(length(ind)>0) {
    x[ind] =  sign(x[ind])*(.Machine$double.xmax/N) 
  }
  
  x0 = mean(.head(x, 2*N), na.rm=TRUE)
  x1 = mean(.tail(x, 2*N), na.rm=TRUE)
  
  dx = x0 - x1
 
  x_tol = reltol * (abs(x0) + reltol)
  
  test = abs(dx) < x_tol
  
  return(test)
}

.smooth_stop3 = function(x, reltol=sqrt(.Machine$double.eps), N=10) {
  
  if(length(x) < 2*N) return(FALSE)
  
  x = .tail(x, 2*N)
  
  ind = which(is.infinite(x))
  if(length(ind)>0) {
    x[ind] =  sign(x[ind])*(.Machine$double.xmax/N) 
  }
  
  x0 = mean(.head(x, N), na.rm=TRUE)
  x1 = mean(.tail(x, N), na.rm=TRUE)
  
  dx = x0 - x1
  
  x_tol = reltol * (abs(x0) + reltol)
  
  test = abs(dx) < x_tol
  
  return(test)
}

.smooth_stop4 = function(x, reltol=sqrt(.Machine$double.eps), N=10) {
  
  if(length(x) < 2*N) return(FALSE)
  
  reltol = N*reltol
  
  x = .tail(x, N)
  
  ind = which(is.infinite(x))
  if(length(ind)>0) {
    x[ind] =  sign(x[ind])*(.Machine$double.xmax/N) 
  }
  
  x0 = max(x, na.rm=TRUE)
  x1 = min(x, na.rm=TRUE)
  
  dx = x0 - x1
  
  x_tol = reltol * (abs(x0) + reltol)
  
  test = dx < x_tol
  
  return(test)
}

.smooth_stop_old = function(x, reltol=sqrt(.Machine$double.eps), N=10, data=FALSE) {
  if(length(x) < 2*N) return(FALSE)
  x = .tail(x, 2*N)
  x[!is.finite(x)] = NA
  dat = data.frame(x=seq_along(x), y=x)
  dat = dat[complete.cases(dat), ]
  n = nrow(dat)
  # mod = scam::scam(y ~ s(x, bs="mpd"), data=dat)
  mod = loess(y ~ x, data=dat)
  dat$value = predict(mod, newdata = dat)
  dat$diff = Inf
  dat$diff[-1] = -diff(dat$value)
  dat$tol = reltol * (abs(dat$value) + reltol)
  if(isTRUE(data)) {
    dat$col = ifelse(dat$diff > dat$tol, "black", "red")
    return(dat)
  }
  return(dat$diff[n] < dat$tol[n])
}

.N_stop = function(x, N) {
  # x = c(na.omit(x))
  if(length(x) <= 2*N) return(FALSE)
  out = mean(.tail(x, N)) > (1-0.1/N)
  return(out)
}

.tail = function(x, n) {
  N = length(x)
  ind = pmax(c(N-n+1, N), 1)
  ind = seq.int(from=ind[1], to=ind[2])
  return(x[ind])
}


.head = function(x, n) {
  N = length(x)
  ind = pmin(c(1, n), N)
  ind = seq.int(from=ind[1], to=ind[2])
  return(x[ind])
}

