.restartCalibration =  function(control, type="restart") {
  res.file = paste0(control$restart.file, ".", type)
  return(file.exists(res.file))
}

.getRestart = function(control, ...) {
  res.file = paste0(control$restart.file, ".restart")
  x = readRDS(res.file)
  msg = sprintf("Restart file '%s' is corrupt.", res.file)
  if(!all(c("opt", "trace") %in% names(x))) stop(msg)
  if(!inherits(x$opt, "ahres.restart")) stop(msg)
  return(x)
}

.getResults = function(control, type="results", ...) {
  res.file = paste0(control$restart.file, ".", type)
  output = readRDS(res.file)
  msg = sprintf("Restart file '%s' is corrupt.", res.file)
  if(!inherits(output, "calibrar.results")) stop(msg)
  class(output) = "list"
  return(output)
}

.createRestartFile = function(opt, trace, control) {
  if(is.null(control$restart.file)) return(invisible())
  if((opt$gen%%control$REPORT)!=0) return(invisible())
  res.file = paste0(control$restart.file, ".restart")
  class(opt) = c("ahres.restart", class(opt))
  force(trace)
  saveRDS(list(opt=opt, trace=trace), file = res.file)
  return(invisible(TRUE))
}

.createRestartFile_Rvmmin = function() {
  return(NULL)
}

.createOutputFile = function(output, control, type="results", phase=0) {
  
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
  
  saveRDS(output, file = res.file)
  
  if(type=="partial") {
    pfile = sprintf(".%s.phase%d", res.file, phase)
    file.copy(from=res.file, to=pfile, overwrite = TRUE)
  }
  
  suppressWarnings(file.remove(paste0(control$restart.file, ".restart"))) 
  return(invisible(NULL))
  
}
