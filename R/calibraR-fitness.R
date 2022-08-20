
#' Calcuted error measure between observed and simulated data
#'
#' @param obs observed data as expected by FUN.
#' @param sim simulated data matching 'obs'
#' @param FUN the error function.
#' @param ... 
#'
#' @return the value of FUN(obs, sim)
#' @export
fitness = function(obs, sim, FUN, ...) {
  FUN = match.fun(FUN)
  output = FUN(obs=obs, sim=sim, ...)
  return(output)
}


# Provided ----------------------------------------------------------------

pois = function(obs, sim, ...) {
  nlogLike = -sum(obs*log(sim) - sim, na.rm=TRUE)
  return(nlogLike)
}

normp = function(obs, sim, ...) {
  penalty = sum((sim)^2, na.rm=TRUE)
  return(penalty)
}

norm2 = function(obs, sim, ...) {
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm2 = function(obs, sim, tiny=1e-2, ...) {
  if(all(!is.finite(sim))) return(Inf)
  obs = log(obs + tiny)
  sim = log(sim + tiny)
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

lnorm3  = function(obs, sim, tiny = 1e-2, ...) {
  if(all(!is.finite(sim))) return(Inf)
  ratio = obs/sim
  ratio[is.nan(ratio)] = NA
  q = mean(ratio, na.rm=TRUE)
  obs = log(obs+tiny) 
  sim = log(sim+tiny)
  nlogLike = sum((obs-sim-log(q))^2, na.rm=TRUE)
  return(nlogLike)
}


multinom = function(sim, obs, size=20, tiny=1e-3) {
  
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
