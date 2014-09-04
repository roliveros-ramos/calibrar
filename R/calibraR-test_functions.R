
# Test functions ----------------------------------------------------------

# method for plotting (2D and 3D)
# method for summary and print
# f() print the minimum

summary.calibrar.function = function(x) {
  "hola"
}


Ackley = function(x) {
  # f(0,0) = 0
  # x,y \in [-5,5]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = -20*exp(-2*sqrt((x^2+y^2)/2)) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y))) + 20 + exp(1)
  return(out) 
}

Sphere = function(x, aggregate=TRUE) {
  # f(0,...,0) = 0
  # x_i \in ]-Inf, Inf[
  out = x^2
  if(isTRUE(aggregate)) return(sum(out)) else return(out) 
}

Rosenbrock = function(x, aggregate=TRUE) {
  # f(1,...,1) = 0
  # x_i \in ]-Inf, Inf[
  out = 100*(x[-1] - x[-length(x)]^2)^2 + (x[-length(x)] - 1)^2
  if(isTRUE(aggregate)) return(sum(out)) else return(out) 
}


Beale = function(x) {
  # f(3, 0.5) = 0
  # x,y \in [-4.5, 4.5]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2
  return(out) 
}


GoldsteinPrice = function(x) {
  # f(0, -1) = 3
  # x,y \in [-2, 2]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  p1 = 19 - 14*x + 3*x^2 -14*y + 6*x*y + 3*y^2
  p2 = 18 - 32*x + 12*x^2 + 48*y - 36*x*y + 27*y^2
  p3 = (x + y + 1)^2
  p4 = (2*x - 3*y)^2
  out = (1 + p3*p1)*(30 + p4*p2)
  return(out) 
}


Booth = function(x) {
  # f(1,3) = 0
  # x,y \in [-10, 10]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = (x + 2*y - 7)^2 + (2*x + y - 5)^2
  return(out) 
}

Bukin6 = function(x) {
  # f(-10, 1) = 0
  # x \in [-15, -5], y \in [-3, 3]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = 100*sqrt(abs(y - 0.01*x^2)) + 0.01*abs(x+10)
  return(out) 
}

Matyas = function(x) {
  # f(0,0) = 0 
  # x,y \in [-10, 10]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = 0.26*(x^2 + y^2) - 0.48*x*y
  return(out) 
}

Levi13 = function(x) {
  # f(1,1) = 0
  # x,y \in [-10, 10]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = sin(3*pi*x)^2 + (x-1)^2*(1 + sin(3*pi*y)^2) + (y-1)^2*(1 + sin(2*pi*y)^2)
  return(out) 
}

ThreeHumpCamel = function(x) {
  # f(0,0) = 0 
  # x,y \in [-5, 5]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = 2*x^2 - 1.05*x^4 + x^6/6 + x*y + y^2
  return(out) 
}

Easom = function(x) {
  # f(pi, pi) = -1
  # x,y \in [-100, 100]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  p1 = (x-pi)^2 + (y-pi)^2
  out = -cos(x)*cos(y)*exp(-p1)
  return(out) 
}

CrossInTray = function(x) {
  # f(+/- 1.34941, +/-1.34941) = -2.06261
  # x,y \in [-10, 10]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = 100 - sqrt(x^2 + y^2)/pi
  out = sin(x)*sin(y)*exp(abs(out))
  out = -0.0001*(abs(out)+1)^0.1
  return(out) 
}

Eggholder = function(x) {
  # f(512, 404.2319) = - 959.6407
  # x,y \in [-512, 512]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = -(y+47)*sin(sqrt(abs(y+0.5*x+47))) - x*sin(sqrt(abs(x-(y+47))))
  return(out) 
}

HolderTable = function(x) {
  # f(+/- 8.05502, +/-9.66459) = -19.2085
  # x,y \in [-10, 10]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = 1 - sqrt(x^2 + y^2)/pi
  out = sin(x)*cos(y)*exp(abs(out))
  out = -abs(out)
  return(out) 
}

McCormick = function(x) {
  # f(-0.54719, -1.54719) = -1.9133
  # x \in [-1.5, 4], y \in [-3, 4]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  out = sin(x+y) + (x-y)^2 - 1.5*x + 2.5*y + 1
  return(out) 
}

Schaffer2 = function(x) {
  # f(0,0) = 0
  # x,y \in [-100, 100]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  p1 = sin(x^2-y^2)^2 - 0.5
  p2 = 1 + 0.001*(x^2 + y^2)
  out = 0.5 + p1/p2^2 
  return(out) 
}

Schaffer4 = function(x) {
  # f(0,1.25313) = 0.292579
  # x,y \in [-100, 100]
  if(length(x)!=2) stop("x must have length 2.")
  y = x[2]
  x = x[1]
  p1 = cos(sin(abs(x^2-y^2))) - 0.5
  p2 = 1 + 0.001*(x^2 + y^2)
  out = 0.5 + p1/p2^2 
  return(out) 
}

StyblinskiTang = function(x, aggregate=TRUE) {
  # f(-2.903534,...,-2.903534) = -3916599n
  # x_i \in [-5, 5]
  out = (x^4 - 16*x^2 + 5*x)/2
  if(isTRUE(aggregate)) return(sum(out)) else return(out) 
}
