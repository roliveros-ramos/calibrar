require(calibrar)
require(optimx)
require(hydroPSO)
require(cmaes)

N = 9 # number of variables in the linear model
T = 100 # number of observations
noise = FALSE # add gaussian noise to the model
shift = FALSE # add a random shift to the slopes
sd = runif(1) # standard deviation of the gaussian noise

# observed data
x = t(matrix(rnorm(N*T, sd=sd), nrow=N, ncol=T))

# slopes for the linear model (real parameters)
slope = seq_len(N) + shift*sample(c(-100, 100), N, replace=TRUE)
# intercept for the linear model (real parameters)
intercept = pi
# real parameters
real = c(intercept, slope)

# function to simulate the linear model
linear = function(x, slope, intercept) {
  stopifnot(length(x)==length(slope))
  out = sum(x*slope) + intercept
  return(out)
}

# simulated data
y = apply(x, 1, linear, slope=slope, intercept=intercept) + noise*rnorm(nrow(x), sd=mean(sd))

# objective function (residual squares sum)
obj = function(par, x, y) {
  intercept = par[1]
  par = par[-1]
  slope = par[seq_len(ncol(x))]
  out = apply(x, 1, linear, slope=slope, intercept=intercept)
  out = sum((out - y)^2)
  return(out)
}

# initial guess for optimization
start = rep(0, N+1)

cat("Running optimization algorithms\n")
cat("\t", date(), "\n")

cat("Running calibrar AHR-ES (unconstrained)\n")
print(system.time(es  <- optimES(par=start, fn=obj, x=x, y=y)))

cat("Running calibrar AHR-ES (constrained)\n")
print(system.time(es2 <- optimES(par=start, fn=obj, x=x, y=y, 
                           lower=rep(-100, length(start)),
                           upper=rep(100, length(start)))))

cat("Running linear model\n")
print(system.time(mod <- lm(y ~ x)))

cat("Running optim CG\n")
print(system.time(opt  <- optim(par=start, fn=obj, x=x, y=y, method="CG")))

cat("Running optim SANN\n")
print(system.time(sann <- optim(par=start, fn=obj, x=x, y=y, method="SANN")))

cat("Running optimx Nelder-Mead + BFGS\n")
print(system.time(optx <- optimx(par=start, fn=obj, x=x, y=y)))

cat("Running hydroPSO\n")
print(system.time(pso  <- hydroPSO(par=start, fn=obj, x=x, y=y, 
               lower=rep(-100, length(start)),
               upper=rep(100, length(start)))))

cat("Running cmaes CMA-ES\n")
print(system.time(cma <- cma_es(par=start, fn=obj, x=x, y=y, 
             lower=rep(-100, length(start)),
             upper=rep(100, length(start)))))


final = rbind(real=real,
              'AHR-ES'=es$par,
              'AHR-ES (constrained)'=es2$par,
              lm=coef(mod),
              SANN=sann$par,
              PSO=pso$par,
              'CMA-ES'=cma$par,
              optx[,seq_along(start)],
              CG=opt$par
              )

print(final)



