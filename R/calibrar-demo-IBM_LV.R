# Generate simulated data -------------------------------------------------

.generateIBMLotkaVolterra = function(path, r=0.09, m=0.045, alpha=5e-4, gamma=0.5,
                                     D=list(N=8e-5, P=8e-5), 
                                     L=list(N=0.3, P=0.3), T=120, replicates=100,
                                     transform=FALSE, seed=880820, ...) {
  
  if(!requireNamespace("ibm", quietly = TRUE)) 
    stop("You need to install the 'ibm' package.")
  
  # 'real' parameters
  
  N0 = as.integer(m/(2*alpha*gamma*L$P))
  P0 = as.integer(r/(2*alpha*L$N))
  
  par_real = list(r=r, m=m, alpha=alpha, gamma=gamma, D=D, L=L, 
                  initial=list(N=N0, P=P0))
  
  set.seed(seed)
  pop = ibm::localLotkaVolterra(par=par_real, T=T, replicates=replicates, 
                                verbose=FALSE, spatial=FALSE, fill=0)
  
  maxpop = .roundmaxpop(max(pop$N, pop$P, na.rm=TRUE), k=5)

  opop = pop
  
  guess = get_guess_LV(pop)
  
  pop$N = rowMeans(pop$N)
  pop$P = rowMeans(pop$P)
  
  main.folder   = file.path(path, "IBMLotkaVolterra")
  data.folder   = file.path(main.folder, "data")
  
  if(!file.exists(data.folder)) dir.create(data.folder, recursive=TRUE)
  
  
  write.csv(pop$N, file.path(data.folder, "N.csv"))
  write.csv(pop$P, file.path(data.folder, "P.csv"))
  
  # parInfo.csv
  
  parInfo = list()
  parInfo$guess = list(r=mean(guess[, "r"], na.rm=TRUE), 
                       m=mean(guess[, "m"], na.rm=TRUE), 
                       alpha=mean(guess[, "alpha"], na.rm=TRUE), 
                       gamma=mean(guess[, "gamma"], na.rm=TRUE), 
                       D=D, L=L,
                       initial=list(N=pop$N[1], P=pop$P[1]))
  
  parInfo$lower = list(r=0, m=0, alpha=1e-8, gamma=1e-8, 
                       D=list(N=1e-8,P=1e-8), L=list(N=1e-4,P=1e-4),
                       initial=list(N=0.5*pop$N[1], P=0.5*pop$P[1]))
  
  parInfo$upper = list(r=1, m=1, alpha=1, gamma=1, 
                       D=list(N=1,P=1), L=list(N=1,P=1),
                       initial=list(N=1.5*pop$N[1], P=1.5*pop$P[1]))
  
  parInfo$phase = list(r=1, m=1, alpha=1, gamma=1, 
                       D=list(N=NA,P=NA), L=list(N=NA,P=NA),
                       initial=list(N=NA, P=NA))
  
  if(isTRUE(transform)) {
   
    xparInfo = parInfo 
    parInfo$guess$alpha = -log10(xparInfo$guess$alpha)
    parInfo$lower$alpha = -log10(xparInfo$upper$alpha)
    parInfo$upper$alpha = -log10(xparInfo$lower$alpha)
    
    parInfo$guess$gamma = .logit(xparInfo$guess$gamma)
    parInfo$lower$gamma = .logit(xparInfo$lower$gamma)
    parInfo$upper$gamma = .logit(xparInfo$upper$gamma)
    
    par_real$alpha     = -log10(par_real$alpha)
    par_real$gamma     = .logit(par_real$gamma)
    
  }
  
  # calibrationInfo.csv
  
  calibrationInfo = list()
  calibrationInfo$variable  = c("N", "P")
  calibrationInfo$type      = "lnorm2"
  calibrationInfo$calibrate = TRUE
  calibrationInfo$weight   = 1
  calibrationInfo$use_data   = TRUE
  calibrationInfo$file   = c("data/N.csv", "data/P.csv")
  calibrationInfo$varid   = "x"
  
  
  calibrationInfo = as.data.frame(calibrationInfo)
  
  write.csv(calibrationInfo, file.path(main.folder, "calibration_settings.csv"), row.names=FALSE)
  
  constants = list(T=T, transform=transform, maxpop=maxpop)
  
  output = c(list(setup=file.path(main.folder, "calibration_settings.csv"), 
                  path = main.folder,
                  par=par_real, bguess=guess, pop=opop), constants, parInfo)
  
  output$value = NA
  
  return(output)
  
}

.roundmaxpop = function(x, k=5) {
  x = 1000*ceiling(k*x/1000)
  x = 10000*ceiling(x/10000)
  return(x)
}


.get_guess_LV = function(N, P, NP, PN) {
  dN = diff(N)
  dP = diff(P)
  
  N2 = head(N, -1) + 0.5*dN
  P2 = head(P, -1) + 0.5*dP
  
  N = head(N, -1) # + 0.5*dN
  P = head(P, -1) # + 0.5*dP
  
  dNr = dN/N
  dPr = dP/P
  
  Px = NP*P # effective number of predators 'predating upon'
  Nx = PN*N # effective number of preys  'being predated'

  P2x = NP*P2 # effective number of predators 'predating upon'
  N2x = PN*N2 # effective number of preys  'being predated'
  
  NP = N*Px
  PN = Nx*P

  NP2 = N2*P2x
  PN2 = N2x*P2
  
  NPg = N*P
  NP2g = N2*P2
  # start point derivative (local interaction)
  modN = lm(dN ~ N + NP + 0)
  modP = lm(dP ~ P + PN + 0)
  
  gamma = pmax(pmin(-coef(modP)[2]/coef(modN)[2], 1), 0)
  
  par = c(coef(modN)[1], -coef(modP)[1], -coef(modN)[2], coef(modP)[2], gamma)
  names(par) = c("r", "m", "alpha", "beta",  "gamma")

  # mid-point derivative (local interaction)
  modN2 = lm(dN ~ N2 + NP2 + 0)
  modP2 = lm(dP ~ P2 + PN2 + 0)
  
  gamma2 = pmax(pmin(-coef(modP2)[2]/coef(modN2)[2], 1), 0)
  
  par0 = c(coef(modN2)[1], -coef(modP2)[1], -coef(modN2)[2], coef(modP2)[2], gamma2)
  names(par0) = c("r0", "m0", "alpha0", "beta0",  "gamma0")

  # start point derivative (global interaction)
  modN3 = lm(dN ~ N + NPg + 0)
  modP3 = lm(dP ~ P + NPg + 0)
  
  gamma3 = pmax(pmin(-coef(modP3)[2]/coef(modN3)[2], 1), 0)
  
  par1 = c(coef(modN3)[1], -coef(modP3)[1], -coef(modN3)[2], coef(modP3)[2], gamma3)
  names(par1) = c("r1", "m1", "alpha1", "beta1",  "gamma1")
  
  # mid- point derivative (global interaction)
  modN4 = lm(dN ~ N2 + NP2g + 0)
  modP4 = lm(dP ~ P2 + NP2g + 0)
  
  gamma4 = pmax(pmin(-coef(modP4)[2]/coef(modN4)[2], 1), 0)
  
  par2 = c(coef(modN4)[1], -coef(modP4)[1], -coef(modN4)[2], coef(modP4)[2], gamma4)
  names(par2) = c("r2", "m2", "alpha2", "beta2",  "gamma2")
  # par0 = c(coef(modN)[1], -coef(modP)[1], -coef(modN)[2], coef(modP)[2], gamma)
  # names(par0) = c("r0", "m0", "alpha0", "beta0",  "gamma0")
  # 
  # modNr = lm(dNr ~ P)
  # modPr = lm(dPr ~ N)
  # gamma = pmax(pmin(-coef(modPr)[2]/coef(modNr)[2], 1), 0)
  # 
  # par1 = c(coef(modNr)[1], -coef(modPr)[1], -coef(modNr)[2], coef(modPr)[2], gamma)
  # names(par1) = c("r1", "m1", "alpha1", "beta1",  "gamma1")
  # 
  # 
  # NPc = NP/coef(modP)[2]
  # modN2 = lm(dN ~ N + NPc + 0)
  # 
  # gamma = -1/coef(modN2)[2]
  # 
  # par2 = c(coef(modN2)[1], -coef(modP)[1], coef(modP)[2]/gamma, coef(modP)[2], gamma)
  # names(par2) = c("r2", "m2", "alpha2", "beta2", "gamma2")
  # 
  par = c(par, par0, par1, par2)
  
  return(par)
  
}

get_guess_LV = function(pop) {
  
  n = ncol(pop$N)
  out = vector("list", n)
  for(i in seq_len(n)) {
    out[[i]] = .get_guess_LV(pop$N[, i], pop$P[, i], pop$NP[, i], pop$PN[, i])
  }
  out = as.data.frame(do.call(rbind, out))
  return(out)
}

.logit = function(p) log(p/(1-p))
