
pois = function(obs, sim, ...) {
  logLike = sum(obs*log(sim) - sim, na.rm=TRUE)
  return(-logLike)
}

normp = function(obs, sim, ...) {
  penalty = sum((sim)^2, na.rm=TRUE)
  return(penalty)
}

norm = function(obs, sim, ...) {
  penalty = sum((obs-sim)^2, na.rm=TRUE)
  return(penalty)
}