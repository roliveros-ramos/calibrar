
commonArgs = function (par, fn, ctrl, rho) {
  rho$n <- n <- length(rho$par <- as.double(par))
  stopifnot(all(is.finite(par)), is.function(fn), length(formals(fn)) >= 
              1)
  rho$.feval. <- integer(1)
  cc <- do.call(function(npt = min(n + 2L, 2L * n), rhobeg = NA, 
                         rhoend = NA, iprint = 0L, maxfun = 10000L, obstop = TRUE, 
                         force.start = FALSE, ...) {
    if (length(list(...)) > 0) 
      warning("unused control arguments ignored")
    list(npt = npt, rhobeg = rhobeg, rhoend = rhoend, iprint = iprint, 
         maxfun = maxfun, obstop = obstop, force.start = force.start)
  }, ctrl)
  ctrl <- new.env(parent = emptyenv())
  lapply(names(cc), function(nm) assign(nm, cc[[nm]], envir = ctrl))
  ctrl$npt <- as.integer(max(n + 2L, min(ctrl$npt, ((n + 1L) * 
                                                      (n + 2L))%/%2L)))
  if (ctrl$npt > (2 * n + 1)) 
    warning("Setting npt > 2 * length(par) + 1 is not recommended.")
  if (is.na(ctrl$rhobeg)) 
    ctrl$rhobeg <- min(0.95, 0.2 * max(abs(par)))
  if (is.na(ctrl$rhoend)) 
    ctrl$rhoend <- 1e-06 * ctrl$rhobeg
  stopifnot(0 < ctrl$rhoend, ctrl$rhoend <= ctrl$rhobeg)
  if (ctrl$maxfun < 10 * n^2) 
    warning("maxfun < 10 * length(par)^2 is not recommended.")
  return(ctrl)
}