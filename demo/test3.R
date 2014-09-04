require(calibrar)

set.seed(12345) # updated to T=20 and L=40 for comparative purposes.

calibrarDemo = function(path=NULL, model=NULL, seed=830613, ...) {
  
  
  
}

.generatePoissonMixedModel = function(path, alpha=0.1, beta=0.05, sd=0.4, T=100, L=12, N0=100, sd_env=0.1, seed=771104, 
                                      ...) {
  
  set.seed(seed)
  
  nsites = ceiling(log10(L+1))
  
  env  = matrix(rnorm(T*L, mean=0, sd=sd_env), nrow=T, ncol=L)
  env  = round(apply(env, 2, cumsum),3)
  
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
  

  parInfo = data.frame(guess = unlist(par_real))
  parInfo$guess = c(0.2, 0.1, rep(0, T-1), sd, log(n[1,]))
  parInfo$min = c(0, 0, rep(-10, T-1), 0.01, rep(0, L))
  parInfo$max = c(1, 0.5, rep(10, T-1), 1, log(1.5*n[1,]))
  parInfo$phase = c(1, 1, rep(2, T-1), NA, rep(3, L))
  
  write.csv(parInfo, file.path(main.folder, "parInfo.csv"))

  # write calibrationInfo.csv
  
  return(par_real)
}

# mode new functions

.PoissonMixedModel = function(par, forcing) {
  # par is a list with 'alpha', 'beta' 'gamma' and 'mu_ini'.
  T = nrow(forcing)
  L = ncol(forcing)
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




