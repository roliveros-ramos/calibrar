
par0 = c(NA, NA)

optimES(par0, fn=Ackley, lower=-5, upper=5)
optimES(par0, fn=Sphere)
optimES(par0, fn=Rosenbrock)
optimES(par0, fn=Beale, lower=-4.5, upper=4.5) #
optimES(par0, fn=GoldsteinPrice, lower=-2, upper=2)
optimES(par0, fn=Booth, lower=-10, upper=10)
optimES(par0, fn=Bukin6, lower=c(-15,-3), upper=c(-5,3))
optimES(par0, fn=Matyas, lower=-10, upper=10)
optimES(par0, fn=Levi13, lower=-10, upper=10)
optimES(par0, fn=ThreeHumpCamel, lower=-5, upper=5)
optimES(par0, fn=Easom, lower=-100, upper=100) #
optimES(par0, fn=CrossInTray, lower=-10, upper=10)
optimES(par0, fn=Eggholder, lower=-512, upper=512) # 
optimES(par0, fn=HolderTable, lower=-10, upper=10)
optimES(par0, fn=McCormick, lower=c(-1.5,-3), upper=4)
optimES(par0, fn=Schaffer2, lower=-100, upper=100)
optimES(par0, fn=Schaffer4, lower=-100, upper=100) #
optimES(par0, fn=StyblinskiTang, lower=-5, upper=5)

# Make it work with multiple phases. 
# 

testx = function(x) Rosenbrock(x) + Sphere(x)
testx2 = function(x) c(Rosenbrock(x), Sphere(x))

optimES(rep(NA,5), fn=Rosenbrock)
