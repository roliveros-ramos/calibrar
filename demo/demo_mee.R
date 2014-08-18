
calibrate(par=rep(NA, 5), fn=Sphere)
calibrate(par=rep(NA, 5), fn=Sphere, replicates=3)
calibrate(par=rep(0.5, 10), fn=Sphere, replicates=3, lower=-5, upper=5)
calibrate(par=rep(0.5, 10), fn=Sphere, replicates=3, lower=-5, upper=5,
          phases=c(1,1,1,2,3))
calibrate(par=rep(0.5, 10), fn=Sphere, replicates=c(1,1,4), lower=-5, upper=5,
          phases=c(1,1,1,2,3))
