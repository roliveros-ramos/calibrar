.restartCalibration =  function(control, type="restart") {
  res.file = paste0(control$restart.file, ".", type)
  return(file.exists(res.file))
}

.getRestart = function(control, ...) {
  res.file = paste0(control$restart.file, ".restart")
  load(res.file)
  if(!exists("opt")) stop("Restart file ", res.file, " is not appropiate.")
  if(!exists("trace")) stop("Restart file ", res.file, " is not appropiate.")
  if(!any(class(opt)=="optimES.restart")) stop("Restart file ", res.file, " is not appropiate.")
  opt   = get("opt")
  trace = get("trace")
  return(list(opt=opt, trace=trace))
}

.getResults = function(control, type="results", ...) {
  res.file = paste0(control$restart.file, ".", type)
  load(res.file)
  if(!exists("output")) stop("Restart file ", res.file, " is not appropiate.")
  if(!any(class(output)=="calibrar.results")) stop("Restart file ", res.file, " is not appropiate.")
  class(output) = "list"
  return(output)
}

.createRestartFile = function(opt, trace, control) {
  if(is.null(control$restart.file)) return(invisible())
  if((opt$gen%%control$REPORT)!=0) return(invisible())
  res.file = paste0(control$restart.file, ".restart")
  class(opt) = c("optimES.restart", class(opt))
  force(trace)
  save(list=c("opt", "trace"), file = res.file, compress=FALSE)
  return(invisible())
}

.createOutputFile = function(output, control, type="results") {
  
  if(is.null(control$restart.file)) return(invisible())
  res.file = paste0(control$restart.file, ".", type)
  class(output) = c("calibrar.results", class(output))
  
  if(type=="results") {
    if(file.exists(res.file)) {
      nfile = sprintf("%s.backup%s", res.file, format(Sys.time(), "%Y%m%d%H%M%S"))
      suppressWarnings(file.rename(from=res.file, to=nfile))
    }
    suppressWarnings(file.remove(paste0(control$restart.file, ".partial")))
  }
  
  save(output, file = res.file, compress=FALSE)
  suppressWarnings(file.remove(paste0(control$restart.file, ".restart"))) 
  return(invisible(NULL))
  
}
