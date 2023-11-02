

.generatePoissonMixedModel = function(path, alpha=0.4, beta=-0.4, T=10, L=6, N0=100, sd_env=0.3, seed=0,
                                      ...) {
  
  nsites = ceiling(log10(L+1))
  
  sd = alpha/2
  
  set.seed(seed)
  env  = matrix(rnorm(T*L, mean=0, sd=sd_env), nrow=T, ncol=L)
  env  = apply(env, 2, cumsum)
  env  = sweep(env, 2, STATS = colMeans(env), FUN = "-")
  env  = env/diff(range(env))/2
  env  = round(env, 4) + 1
  
  colnames(env) = paste0("L", sprintf(paste0("%0", nsites, "d"), seq_len(L)))
  rownames(env) = seq_len(T)
  
  # 'real' parameters
  gamma  = rnorm(T-1, mean=0, sd=sd)
  mu_ini = log(rpois(n=L, lambda=runif(L, min=N0/5, max=2*N0))) #initial means
  
  par_real = list(alpha=alpha, beta=beta, gamma=gamma, sd=sd, mu_ini=mu_ini)
  
  mu = .PoissonMixedModel(par=par_real, forcing=env)
  mu = as.matrix(as.data.frame(mu))
  # observed abundances
  n = matrix(rpois(length(mu), lambda=mu), nrow=T, ncol=L)  
  
  main.folder   = file.path(path, "PoissonDemo")
  data.folder   = file.path(main.folder, "data")
  master.folder = file.path(main.folder, "master")
  
  if(!file.exists(data.folder)) dir.create(data.folder, recursive=TRUE)
  if(!file.exists(master.folder)) dir.create(master.folder, recursive=TRUE)
  
  write.csv(env, file.path(master.folder, "environment.csv"))
  
  for(i in seq_len(L)) {
    isite = sprintf(paste0("%0", nsites, "d"), i)
    site.file = paste0("site_", isite, ".csv")
    dat = matrix(n[, i], ncol=1)
    colnames(dat) = paste0("L",isite)
    write.csv(dat, file.path(data.folder, site.file))
  }
  
  # parInfo.csv
  
  parInfo = list()
  parInfo$guess = relist(c(0.2, 0.1, rep(0, T-1), par_real$sd, round(log(n[1,]),3)), par_real)
  parInfo$lower = relist(c(0, -0.5, rep(-5*sd, T-1), 0, round(rep(min(log(0.5*n[1,])),L),1)), par_real)
  parInfo$upper = relist(c(1, 0.5, rep(5*sd, T-1), 1, round(rep(max(log(1.5*n[1,])),L),1)), par_real)
  parInfo$phase = relist(c(1, 1, rep(1, T-1), NA, rep(1, L)), par_real)
  
  # calibrationInfo.csv
  
  calibrationInfo = list()
  calibrationInfo$variable   = paste0("site_", sprintf(paste0("%0", nsites, "d"), seq_len(L)))
  calibrationInfo$type       = "pois"
  calibrationInfo$calibrate  = TRUE
  calibrationInfo$weight     = 1
  calibrationInfo$use_data   = TRUE
  calibrationInfo$file       = file.path("data", paste0(calibrationInfo$variable, ".csv"))
  calibrationInfo$varid      = paste0("L", sprintf(paste0("%0", nsites, "d"), seq_len(L)))
  
  calibrationInfo = as.data.frame(calibrationInfo)
  
  calibrationInfo2 = list()
  calibrationInfo2$variable  = "gammas"
  calibrationInfo2$type      = "normp"
  calibrationInfo2$calibrate = TRUE
  calibrationInfo2$weight    = 1/(2*sd^2)
  calibrationInfo2$use_data  = FALSE
  calibrationInfo2$file      = NA 
  calibrationInfo2$varid     = NA
  
  calibrationInfo2 = as.data.frame(calibrationInfo2)
  calibrationInfo = rbind(calibrationInfo, calibrationInfo2)
  
  write.csv(calibrationInfo, file.path(main.folder, "calibration_settings.csv"),
            row.names=FALSE)
  
  constants = list(T=T, L=L)
  
  output = c(list(path=main.folder, par=par_real, 
                  setup=file.path(main.folder, "calibration_settings.csv")), constants, parInfo)
 
  setup = calibration_setup(file = output$setup)
  observed = calibration_data(setup=setup, path=output$path)
  forcing = as.matrix(read.csv(file.path(output$path, "master", "environment.csv"), row.names=1))
  
  run_model = function(par, forcing) {
    output = calibrar:::.PoissonMixedModel(par=par, forcing=forcing)
    output = c(output, list(gammas=par$gamma)) # adding gamma parameters for penalties
    return(output)
  }
  
  obj2 = calibration_objFn(model=run_model, setup=setup, observed=observed, 
                           forcing=forcing, aggregate=TRUE)
  
  output$value = obj2(output$par)
  
  return(output)
  
}

# mode new functions

.PoissonMixedModel = function(par, forcing) {
  # par is a list with 'alpha', 'beta' 'gamma', 'sd' and 'mu_ini'.
  T = nrow(forcing)
  L = ncol(forcing)
  alpha  = par$alpha
  beta   = par$beta
  gamma  = if(!is.null((par$gamma))) par$gamma else rep(0, T-1)
  # gamma  = gamma - mean(gamma)
  
  mu = matrix(exp(par$mu_ini), nrow=T, ncol=L, byrow = TRUE)
  
  # mu[1,] = mu_ini
  
  rate = exp(apply(rbind(0, alpha + beta*forcing[-T, ] + gamma), 2, cumsum))
  mu   = mu*rate 

  # for(t in seq_len(T-1)) {
  #   # log_mu_new = log(mu[t,]) + alpha + beta*forcing[t,] + gamma[t]
  #   # log_mu_new = log(mu[t,]) + r0[t,]
  #   # mu[t+1, ] = exp(log_mu_new)
  #   mu[t+1, ] = mu[t, ]*r0[t, ]
  # }
  
  output = as.list(as.data.frame(mu)) # return a list with the results
  # names of the outputs matching observed data names
  names(output) = paste0("site_", sprintf(paste0("%0", ceiling(log10(L+1)), "d"), seq_len(L)))
  
  return(output)
}
