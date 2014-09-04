
calibrate(par=rep(NA, 5), fn=Sphere)
calibrate(par=rep(NA, 5), fn=Sphere, replicates=3)
calibrate(par=rep(0.5, 10), fn=Sphere, replicates=3, lower=-5, upper=5)
calibrate(par=rep(0.5, 10), fn=Sphere, replicates=3, lower=-5, upper=5,
          phases=c(1,1,1,2,3))
calibrate(par=rep(0.5, 10), fn=Sphere, replicates=c(1,1,4), lower=-5, upper=5,
          phases=c(1,1,1,2,3))


# Appendix 2

# data/variable_n.csv
# master/.
# parInfo.csv (opt, min, max, phase)
# calibrationInfo.csv (...)

# model1_run/i0
# model1_restart.rda # restart file for model1
# model1_outputs.csv

# TRACE FOR LATER!

# Create demostration data
runDemo(path=NULL)



runModel = function(par, ...) {
  
  # organizing parameters
  K1 = 1 # from par
  K2 = 1
  r1 = 0.3
  r2 = 0.2
  h  = 0.5
  N0 = 0.2
  output = .model(N0=N0, K=c(K1, K2), r=c(r1,r2), h=h)
  
}


