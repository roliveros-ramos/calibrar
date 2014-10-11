.restartCalibration =  function(control) {
  res.file = paste0(control$restart.file, ".restart")
  return(file.exists(res.file))
}

.getRestart = function(control=control, ...) {
  res.file = paste0(control$restart.file, ".restart")
  load(res.file)
  if(!exists(opt)) stop("Restart file is not appropiate.")
  if(!exists(trace)) stop("Restart file is not appropiate.")
  if(class(opt)!="optimES.restart") stop("Restart file is not appropiate.")
  return(list(opt=opt, trace=trace))
}


.createRestartFile = function(opt, trace=trace, control) {
  res.file = paste0(control$restart.file, ".restart")
  class(opt) = c("optimES.restart", class(opt))
  save(list=c("opt", "trace"), file = res.file)
  return(invisible())
}
