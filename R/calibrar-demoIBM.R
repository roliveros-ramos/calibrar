# Generate simulated data -------------------------------------------------

.generateIBMLotkaVolterra = function(path, r=0.1, m=0.05, alpha=5e-4, beta=5e-4,
                                     D=list(N=8e-5, P=8e-5), 
                                     L=list(N=0.2, P=0.2), T=100, 
                                     uselog=FALSE, seed=880820, ...) {
  
  if(!requireNamespace("ibm", quietly = TRUE)) 
    stop("You need to install the 'ibm' package.")
  
  # 'real' parameters
  
  N0 = m/(2*beta*L$P)
  P0 = r/(2*alpha*L$N)
  
  par_real = list(r=r, m=m, alpha=alpha, beta=beta, D=D, L=L, 
                  initial=list(N=N0, P=P0))
  
  set.seed(seed)
  pop = ibm:::.localLotkaVolterra(par=par_real, T=T, verbose=FALSE, spatial=FALSE)
  
  main.folder   = file.path(path, "IBMLotkaVolterra")
  data.folder   = file.path(main.folder, "data")
  
  if(!file.exists(data.folder)) dir.create(data.folder, recursive=TRUE)
  
  write.csv(pop$N, file.path(data.folder, "N.csv"))
  write.csv(pop$P, file.path(data.folder, "P.csv"))
  
  # parInfo.csv
  
  parInfo = list()
  parInfo$guess = list(r=0.1, m=0.1, alpha=1e-3, beta=1e-3, 
                       D=list(N=8e-5, P=8e-5), L=list(N=0.2, P=0.2),
                       initial=list(N=pop$N[1], P=pop$P[1]))
  parInfo$lower = list(r=0, m=0, alpha=1e-8, beta=1e-8, 
                       D=list(N=1e-8,P=1e-8), L=list(N=1e-4,P=1e-4),
                       initial=list(N=0.5*pop$N[1], P=0.5*pop$P[1]))
  parInfo$upper = list(r=1, m=1, alpha=1, beta=1, 
                       D=list(N=1,P=1), L=list(N=1,P=1),
                       initial=list(N=1.5*pop$N[1], P=1.5*pop$P[1]))
  
  parInfo$phase = list(r=1, m=1, alpha=1, beta=1, 
                       D=list(N=NA,P=NA), L=list(N=NA,P=NA),
                       initial=list(N=NA, P=NA))
  
  if(isTRUE(uselog)) {
   
    xparInfo = parInfo 
    parInfo$guess$alpha = -log10(xparInfo$guess$alpha)
    parInfo$lower$alpha = -log10(xparInfo$upper$alpha)
    parInfo$upper$alpha = -log10(xparInfo$lower$alpha)
    
    parInfo$guess$beta = -log10(xparInfo$guess$beta)
    parInfo$lower$beta = -log10(xparInfo$upper$beta)
    parInfo$upper$beta = -log10(xparInfo$lower$beta)
    
    par_real$alpha     = -log10(par_real$alpha)
    par_real$beta      = -log10(par_real$beta)
    
  }
  
  # calibrationInfo.csv
  
  calibrationInfo = list()
  calibrationInfo$variable  = c("N", "P")
  calibrationInfo$type      = "lnorm2"
  calibrationInfo$calibrate = TRUE
  calibrationInfo$weights   = 1
  calibrationInfo$useData   = TRUE
  
  calibrationInfo = as.data.frame(calibrationInfo)
  
  write.csv(calibrationInfo, file.path(main.folder, "calibrationInfo.csv"), row.names=FALSE)
  
  constants = list(T=T, uselog=uselog)
  
  output = c(list(path=main.folder, par=par_real), constants, parInfo)
  
  return(output)
  
}


