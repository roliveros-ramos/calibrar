.generateSIRModel = function(path, beta=0.4, gamma=0.2, T=100, 
                                      S0=999, I0=1, R0=0, seed=0, ...) {

  set.seed(seed)
  # 'real' parameters
  par_real = list(beta=beta, gamma=gamma, initial=list(S=S0, I=I0, R=R0))
  
  n = .SIRModel(par=par_real, T=T)
  
  # observed abundances
  # n = rapply(pop, f=jitter, how = "list")
  
  main.folder   = file.path(path, "SIRDemo")
  data.folder   = file.path(main.folder, "data")
  
  if(!file.exists(data.folder)) dir.create(data.folder, recursive=TRUE)
    
  for(i in c("susceptible", "infected", "recovered")) {
    ifile = paste0(i, ".csv")
    dat = matrix(n[[i]], ncol=1)
    colnames(dat) = i
    write.csv(dat, file.path(data.folder, ifile))
  }
  
  # parInfo.csv
  ini = c(n$susceptible[1], n$infected[1], n$recovered[1])
  lower = pmax(floor(0.5*(ini-1)), 0)
  upper = pmax(ceiling(1.5*(ini+1)), 0)
  
  parInfo = list()
  parInfo$guess = list(beta=0.1, gamma=0.1, initial=list(S=ini[1], I=ini[2], R=ini[3]))
  parInfo$lower = list(beta=0, gamma=0, initial=list(S=lower[1], I=lower[2], R=lower[3]))
  parInfo$upper = list(beta=100, gamma=100, initial=list(S=upper[1], I=upper[2], R=upper[3]))
  parInfo$phase = list(beta=1, gamma=1, initial=list(S=NA, I=NA, R=NA))
  
  # calibrationInfo.csv
  
  calibrationInfo = list()
  calibrationInfo$variable  = c("susceptible", "infected", "recovered")
  calibrationInfo$type      = "lnorm2"
  calibrationInfo$calibrate = TRUE
  calibrationInfo$weight   = 1
  calibrationInfo$use_data   = TRUE
  calibrationInfo$file   = c("data/susceptible.csv", "data/infected.csv", "data/recovered.csv")
  calibrationInfo$varid   = c("susceptible", "infected", "recovered")
  
  calibrationInfo = as.data.frame(calibrationInfo)
  
  write.csv(calibrationInfo, file.path(main.folder, "calibration_settings.csv"), row.names=FALSE)

  constants = list(T=T)
  
  output = c(list(setup=file.path(main.folder, "calibration_settings.csv"), 
                  path = main.folder,
                  par=par_real), constants, parInfo)
  
  setup = calibration_setup(file = output$setup)
  observed = calibration_data(setup=setup, path=output$path)
  
  run_model = .SIRModel
  
  obj2 = calibration_objFn(model=run_model, setup=setup, observed=observed, T=T, 
                           aggregate=TRUE)
  
  output$value = obj2(output$par)
  
  return(output)
  
}



.SIRModel = function(par, T) {
  if(!requireNamespace("deSolve", quietly = TRUE))
    stop("You need to install the 'deSolve' package.")
  # par is a list with 'alpha', 'beta' 'gamma', 'sd' and 'mu_ini'.
  SIR = function(t, y, parms, ...) {
    N = sum(unlist(parms$initial))
    beta = parms$beta
    gamma = parms$gamma
    S = y[1]
    I = y[2]
    dS = -beta*S*I/N
    dI = +beta*S*I/N -gamma*I
    dR = +gamma*I
    return(list(c(dS, dI, dR)))
  }
  times = seq(0, T)
  y0 = c(par$initial$S, par$initial$I, par$initial$R)
  sol = deSolve::ode(y=y0, times=times, func=SIR, parms=par, method="ode45")
  out = as.list(as.data.frame(sol[,-1]))
  names(out) = c("susceptible", "infected", "recovered")
  out$susceptible[is.na(out$susceptible)] = 0
  out$infected[is.na(out$infected)] = 0
  out$recovered[is.na(out$recovered)] = 0
  return(out)
}
