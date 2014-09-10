

.generatePoissonMixedModel = function(path, alpha=0.4, beta=0.4, T=10, L=6, N0=100, sd_env=0.3, 
                                      ...) {
  
  nsites = ceiling(log10(L+1))
  
  sd = alpha
  
  env  = matrix(rnorm(T*L, mean=0, sd=sd_env), nrow=T, ncol=L)
  env  = env/diff(range(env))
#   env  = apply(env, 2, filter, filter=rep(1,3)/3)
#   env  = env/sqrt(seq_len(T))
  env  = round(env, 4)
  
  colnames(env) = paste0("L", sprintf(paste0("%0", nsites, "d"), seq_len(L)))
  rownames(env) = seq_len(T)
  
  # 'real' parameters
  gamma  = rnorm(T-1, mean=0, sd=sd)
  mu_ini = log(rpois(n=L, lambda=runif(L, min=N0/4, max=N0))) #initial means
  
  par_real = list(alpha=alpha, beta=beta, gamma=gamma,  
                  sd=sd, mu_ini=mu_ini)
  
  mu = .PoissonMixedModel(par=par_real, forcing=env)
  
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
  
  parInfo = data.frame(guess = unlist(par_real))
  parInfo$guess = c(0.2, 0.1, rep(0, T-1), alpha/3, round(log(n[1,]),3))
  parInfo$lower = c(0, -0.5, rep(-1.5, T-1), 0, rep(0, L))
  parInfo$upper = c(1, 0.5, rep(1.5, T-1), 1, round(rep(max(log(1.5*n[1,])),L),1))
  parInfo$phase = c(1, 1, rep(2, T-1), NA, rep(3, L))
  
  write.csv(parInfo, file.path(main.folder, "parInfo.csv"))
  
  # calibrationInfo.csv
  
  calibrationInfo = list()
  calibrationInfo$variable  = paste0("site_", sprintf(paste0("%0", nsites, "d"), seq_len(L)))
  calibrationInfo$type      = "pois"
  calibrationInfo$calibrate = TRUE
  calibrationInfo$weights   = 1
  calibrationInfo$useData   = TRUE
  
  calibrationInfo = as.data.frame(calibrationInfo)
  
  calibrationInfo2 = list()
  calibrationInfo2$variable = "gammas"
  calibrationInfo2$type     = "normp"
  calibrationInfo2$calibrate = TRUE
  calibrationInfo2$weights   = 1/(2*sd^2)
  calibrationInfo2$useData   = FALSE
  
  calibrationInfo2 = as.data.frame(calibrationInfo2)
  calibrationInfo = rbind(calibrationInfo, calibrationInfo2)
  
  write.csv(calibrationInfo, file.path(main.folder, "calibrationInfo.csv"),
            row.names=FALSE)
  
  constants = list(T=T, L=L)
  
  output = list(path=main.folder, par=par_real, constants)
  
  return(output)
  
}

# mode new functions

.PoissonMixedModel = function(par, forcing) {
  # par is a list with 'alpha', 'beta' 'gamma', 'sd' and 'mu_ini'.
  T = nrow(forcing)
  L = ncol(forcing)
  forcing = as.matrix(forcing)
  alpha  = par$alpha
  beta   = par$beta
  gamma  = if(!is.null((par$gamma))) par$gamma else rep(0, T-1)
  mu_ini = exp(par$mu_ini)
  
  mu = matrix(nrow=T, ncol=L)
  
  mu[1,] = mu_ini
  
  for(t in seq_len(T-1)) {
    log_mu_new = log(mu[t,]) + alpha + beta*forcing[t,] + gamma[t]
    mu[t+1, ] = exp(log_mu_new)
  }
  return(mu)
}

