
# wrapper for cma_es ------------------------------------------------------

.cmaes = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  npar = length(unlist(par))
  
  if(!is.null(control$trace)) control$trace = ifelse(control$trace>1, TRUE, FALSE)
  
  # all default values!
  con = list(trace = FALSE, fnscale = 1, stopfitness = -Inf, maxit = 100 * npar^2,
             sigma = 0.5, keep.best = TRUE, vectorized = FALSE, diag = FALSE,
             lambda = 4 + floor(3 * log(npar)))
  
  con = within(con, {
    stop.tolx = 1e-12 * sigma
    diag.sigma = diag
    diag.eigen = diag
    diag.value = diag
    diag.pop = diag
    mu = floor(lambda/2)
    weights = log(mu + 1) - log(1:mu)
    weights = weights/sum(weights)
    mueff = sum(weights)^2/sum(weights^2)
    ccum = 4/(npar + 4)
    cs = (mueff + 2)/(npar + mueff + 3)
    ccov.mu = mueff
    ccov.1 = (1/ccov.mu) * 2/(npar + 1.4)^2 + (1 - 1/ccov.mu) * ((2 * ccov.mu - 1)/((npar + 2)^2 + 2 * ccov.mu))
    damps = 1 + 2 * max(0, sqrt((mueff - 1)/(npar + 1)) - 1) + cs
  })
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(cmaes::cma_es(par=par, fn=fn, lower=lower, upper=upper, control=control))
  
  if(is.null(output$par)) {
    output$par = relist(rep(NA, npar), skeleton=par)
    if(!is.finite(output$value)) warning("Infinite value reached for fn.")
  }
  
  return(output)
  
}



# wrapper for genSA -------------------------------------------------------

.genSA = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  con = list(maxit = 5000, threshold.stop = NULL, temperature = 5230, 
             visiting.param = 2.62, acceptance.param = -5, max.time = NULL, 
             nb.stop.improvement = 1e+06, smooth = TRUE, max.call = 1e+07, 
             simple.function = FALSE, trace.fn = NULL, verbose = FALSE, 
             trace.mat = TRUE, seed = -100377, high.dim = TRUE, tem.restart = 0.1,
             markov.length = 2 * length(lower))
  
  control = check_control(control=control, default=con)
  
  xoutput = suppressWarnings(GenSA::GenSA(par=par, fn=fn, lower=lower, 
                                         upper=upper, control=control))
  
  output = list()
  output$par  = xoutput$par
  output$value = xoutput$value
  output$counts = c('function'=xoutput$counts, 'gradient'=0)	
  output$convergence = 0
  output$message = NA
  output$hessian = NULL
  
  return(output)
  
}



# wrapper for DEoptim -----------------------------------------------------

.DE = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(isTRUE(control$parallel)) control$parallelType="auto"
  
  con = DEoptim::DEoptim.control()
  if(is.null(control$NP)) control$NP = 10*length(lower)
  if(is.na(control$NP)) control$NP = 10*length(lower)
  if(control$NP < 4) control$NP = 10*length(lower)
  if(!is.null(control$trace))  control$trace = ifelse(control$trace>=1, TRUE, FALSE)
  if(is.null(control$trace)) control$trace = FALSE
  
  control = check_control(control=control, default=con)
  
  # we pass 'par' as initial population (gaussian with uniform-like variance).
  if(is.null(control$initialpop)) {
    control$initialpop = t(.createRandomPopulation(n=control$NP, mean=par, lower=lower, upper=upper)) 
  }

  # think how to pass 'par'. Initial pop?
  xoutput = suppressWarnings(DEoptim::DEoptim(fn=fn, lower=lower, upper=upper, control=control))
  
  output = list()
  output$par  = xoutput$optim$bestmem
  output$value = xoutput$optim$bestval
  output$counts = c('function'=xoutput$optim$nfeval,
                    gradient=NA)
  output$convergence = NA
  output$message = NULL
  
  return(output)
  
}
# wrapper for soma --------------------------------------------------------

.soma = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  control = NULL # check!
  
  con = soma::all2one()
  if(is.null(control$populationSize)) control$populationSize = con$populationSize
  
  reportr::setOutputLevel(6)
  if(isTRUE(control$verbose)) reportr::setOutputLevel(3) 
  control$verbose = NULL
  
  control = check_control(control=control, default=con)

  # we pass 'par' as initial population (gaussian with uniform-like variance).
  init = .createRandomPopulation(n=control$populationSize, mean=par, lower=lower, upper=upper) 
  
  xoutput = suppressWarnings(soma::soma(costFunction=fn, 
                                        bounds=list(min=lower, max=upper), 
                                        options=control, init=init))
  
  output = list()
  output$par  = xoutput$population[, xoutput$leader]
  output$value = xoutput$cost[xoutput$leader]
  output$counts = c('function'=xoutput$migrations*length(xoutput$cost),
                    gradient=NA)
  output$convergence = NA
  output$message = NULL
  
  return(output)
  
}


# wrapper for genoud ------------------------------------------------------

.genoud = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(is.null(control$print.level))
    control$print.level = if(isTRUE(control$verbose)) control$trace else 0
  if(!is.null(control$maxit)) control$max.generations = control$maxit
    
  con = list(pop.size=1000, max.generations=100, 
             wait.generations=10, hard.generation.limit=TRUE, 
             MemoryMatrix=TRUE, default.domains=10, 
             solution.tolerance=0.001, boundary.enforcement=0, lexical=FALSE,
             gradient.check=TRUE, BFGS=TRUE, data.type.int=FALSE,  
             unif.seed=round(runif(1, 1, 2147483647L)), 
             int.seed=round(runif(1, 1, 2147483647L)), share.type=0, print.level=2,
             instance.number=0, output.path="stdout", output.append=FALSE, 
             project.path=NULL, P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, 
             P8=50, P9=0, P9mix=NULL, BFGSburnin=0, BFGSfn=NULL, BFGShelp=NULL, 
             control=list(), transform=FALSE, debug=FALSE, cluster=FALSE, balance=FALSE)
  con$optim.method=ifelse(con$boundary.enforcement < 2, "BFGS", "L-BFGS-B")
  
  control = check_control(control=control, default=con)
  
  control$starting.values = par
  control$fn              = fn
  control$gr              = gr
  control$nvars           = length(par) 
  control$max             = FALSE 
  control$Domains         = cbind(lower, upper)
  control$hessian         = hessian 
  
  xoutput = suppressWarnings(do.call(rgenoud::genoud, args=control))
  
  output = list()
  output$par  = xoutput$par
  output$value = xoutput$value
  output$counts = c('function'=xoutput$popsize*xoutput$generations, gradient=NA)
  output$convergence = NA
  output$message = NULL
  output$hessian = xoutput$hessian
  
  return(output)
  
}

# wrapper for pso ---------------------------------------------------------

.pso = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  hybrid = if(method=="hybridPSO") TRUE else FALSE
  
  if(!is.null(control$type)) control$pso.method = control$type
    
  method = if(!is.null(control$pso.method)) control$pso.method else "SPSO2007"
  # opso = c("PSO", "PSO2007", "SPSO2007", "hybridPSO")
  # method = if(method %in% opso) "SPSO2007" else "SPSO2011"
  control$type   = method 
  if(hybrid) control$hybrid = hybrid
  
  con = list(trace = 0, maxit = 1000, maxf = Inf, abstol = -Inf, reltol = 0, 
             REPORT = 10L, trace.stats = FALSE, 
             s = if(method=="SPSO2011") 40 else floor(10+2*sqrt(length(par))),
             k = 3, w = 1/(2*log(2)), c.p = 0.5+log(2), c.g = 0.5+log(2),
             d = NULL, v.max = NA, rand.order = TRUE, max.restart = Inf,
             maxit.stagnate = Inf, vectorize = FALSE, hybrid = FALSE, 
             hybrid.control = NULL, type = "SPSO2007")
  con$p = 1-(1-1/con$s)^con$k
  
  control = check_control(control=control, default=con)
  
  output = suppressWarnings(pso::psoptim(par=par, fn=fn, gr=gr, lower=lower, 
                                         upper=upper, control=control))
  
  return(output)
  
}




