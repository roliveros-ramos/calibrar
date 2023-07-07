# all methods

.all_methods = function(method, par, fn, gr, lower, upper, control, hessian) {
  output = 
    switch(method, 
           "Nelder-Mead" = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="Nelder-Mead"), 
           "BFGS"        = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="BFGS"), 
           "CG"          = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="CG"), 
           "L-BFGS-B"    = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="L-BFGS-B"),
           "SANN"        = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="SANN"),
           "Brent"       = .optim(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian, method="Brent"), 
           "nlm"         = .nlm(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "nlminb"      = .nlminb(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "Rcgmin"      = .Rcgmin(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "Rvmmin"      = .Rvmmin(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "hjn"         = .hjn(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "spg"         = .spg(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "LBFGSB3"     = .lbfgsb3(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           "AHR-ES"      = .ahres(par=par, fn=fn, gr=gr, lower=lower, upper=upper, control=control, hessian=hessian), 
           stop(printf("UNSUPPORTED METHOD: %s.", sQuote(method)), call. = FALSE)
    )
  return(output)
}

# Manage control options --------------------------------------------------

check_control = function(control, default, minimal=TRUE, verbose=TRUE) {
  
  control = control[!sapply(control, is.null)] # remove NULL values
  nm_full = names(default)
  ignored = setdiff(names(control), nm_full)
  keep = names(default)[nm_full %in% names(control)]
  if(isTRUE(minimal)) return(control[keep])
  msg = sprintf("Ignoring control arguments: %s.", paste(sQuote(ignored), collapse=", "))
  if(length(ignored) & verbose) message(msg)
  default[keep] = control[keep]
  return(default)
}


# Auxiliar functions to run fn on disk ------------------------------------

copy_master_folder = function(control, n=NULL) {
  if(is.null(control$master)) return(invisible())
  if(is.null(control$run))
    stop("You must specify a 'run' directory to copy the contents of the 'master' folder.")
  if(is.null(n)) return(invisible())
  for(i in (seq_len(n+1) - 1)) .copyMaster(control, i)
  return(invisible())
}

.copyMaster = function(control, i) {
  newDir = .getWorkDir(control$run, i)
  if(!dir.exists(newDir)) dir.create(newDir, recursive=TRUE)
  if(isTRUE(control$master)) return(invisible(newDir))
  if(!dir.exists(control$master)) stop("The 'master' directory does not exist.") 
  xfiles = dir(path=control$master)
  from   = file.path(control$master, xfiles)
  file.copy(from=from, to=newDir, recursive=TRUE, overwrite = FALSE) # why?
  return(invisible(newDir))
}

.getWorkDir = function(run, i) {
  if(is.null(run)) return(invisible())
  work.dir = file.path(run, paste0("i", i))
  return(work.dir)
}
