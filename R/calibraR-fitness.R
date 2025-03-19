
#' Calcuted error measure between observed and simulated data
#'
#' @param obs observed data as expected by FUN.
#' @param sim simulated data matching 'obs'
#' @param FUN the error function. Current accepted values area: 'norm2', 
#' 'lnorm2', 'lnorm3', 'multinomial', 'pois', 'penalty0', 'penalty1', 'penalty2' and 'normp'.
#' @param ... Additional arguments to FUN
#'
#' @return the value of FUN(obs, sim, ...)
#' @export
objFn = function(obs, sim, FUN, ...) {
  FUN = match.fun(FUN)
  output = FUN(obs=obs, sim=sim, ...)
  return(output)
}

#' @export
#' @rdname objFn 
fitness = objFn

# Penalties ---------------------------------------------------------------

normp = function(obs, sim, ...) {
  penalty = sum((sim)^2, na.rm=TRUE)
  return(penalty)
}

re = normp

penalty = function(obs, sim, n=100, ...) {
  # assumes a fixed sample size of 'n'
  penalty = n*mean((sim)^2, na.rm=TRUE)
  return(penalty)
}

# identical to penalty, used for testing.
penalty2 = penalty


# Likelihoods -------------------------------------------------------------

pois = function(obs, sim, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  nlogLike = -sum(obs*log(sim) - sim, na.rm=TRUE)
  return(nlogLike)
}

norm2 = function(obs, sim, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm2 = function(obs, sim, tiny=1e-2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  obs = log(obs + tiny)
  sim = log(sim + tiny)
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm3  = function(obs, sim, tiny = 1e-2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  ratio = obs/sim
  ratio[is.nan(ratio)] = NA
  q = mean(ratio, na.rm=TRUE)
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm4  = function(obs, sim, tiny = 1e-2, b=1, c=2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  penq = rangeq(obs=obs, sim=sim, b=b, c=c, dump=TRUE)
  q = attr(penq, "q")
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike + penq)
}

lnorm4b  = function(obs, sim, tiny = 1e-2, b=1, c=2, ...) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  penq = rangeq(obs=obs, sim=sim, b=b, c=c, dump=FALSE)
  q = attr(penq, "q")
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike + penq)
}

rangeq = function(obs, sim, b=1, c=2, dump=TRUE) {
  if(all(is.na(obs))) stop("All observed values are NA.")
  ratio = obs/sim
  ratio[is.nan(ratio)] = NA
  n = sum(!is.na(ratio))
  q = mean(ratio, na.rm=TRUE)
  if(isTRUE(dump)) {
    pen = n*(pmax(abs((log2(q))), b)^c - b^c)
  } else {
    pen = sum((pmax(abs((log2(ratio))), b)^c - b^c))
  }
  attr(pen, "q") = q
  return(pen)
}

multinom = function(sim, obs, size=20, tiny=1e-3) {
  
  if(all(is.na(obs))) stop("All observed values are NA.")
  
  A = ncol(sim) # number of classes
  
  # checking for simulated values with only zeros and replace ones (assumed uniform)
  sim.allzero = apply(sim, 1, FUN=function(x) all(na.omit(x) == 0) & !all(is.na(x)) )
  sim[which(sim.allzero), ] = sim[which(sim.allzero), ] + 1
  
  sim.sum = rowSums(sim, na.rm=TRUE) # only zero for all NAs.
  obs.sum = rowSums(obs, na.rm=TRUE)
  
  # setting to NA for numerical convenience
  sim.sum[sim.sum==0] = NA
  
  # removing 'all zeros' from obs, because it means we have no proportion data
  obs.sum[which(obs.sum==0)] = NA  
  
  Psim     = sim/sim.sum
  Pobs     = obs/obs.sum
  
  sigma2 = ((1-Pobs)*Pobs + 1/A)/size
  
  error = log(exp(-((Pobs - Psim)^2)/(2*sigma2)) + tiny)
  
  nlogLike = -size*sum(error, na.rm=TRUE)
  return(nlogLike)
}
