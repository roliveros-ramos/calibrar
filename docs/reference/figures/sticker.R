
# RStudio: https://github.com/rstudio/hex-stickers
# hexbin: http://hexb.in/
# hexbin github: https://github.com/maxogden/hexbin
# make stickers with: https://github.com/GuangchuangYu/hexSticker


# Note the hexbin github readme says to use
# 181x209 as a png (preview on mac can easily resize)
# and also provide an svg
# I used this to convert png to svg: https://www.aconvert.com/image/png-to-svg/


# https://coolors.co/193763-f4f4f4-00173d-ff5a5f-c81d25

library(hexSticker)
library(ggplot2)
library(gts)
library(calibrar)

mt = expand.grid(x=seq_len(ncol(volcano)), 
                 y=seq_len(nrow(volcano)))
mt$z = as.numeric(t(volcano))
mt$x = (mt$x-1.01)/(max(mt$x)-1)
mt$y = (mt$y-1.01)/(max(mt$y)-1)
xout = seq(0.01, 0.99,length.out=1000)

z = interpolate(mt$x, mt$y, mt$z, xout=xout, yout=xout)
# z$z = t(z$z)

mt2 = expand.grid(x=xout, y=xout)
mt2$z = as.numeric(z$z)

z$xb = gts:::.getBreaks(z$x)
z$yb = gts:::.getBreaks(z$y)

find_top = function(par, z) {
  x = par[1]
  y = par[2]
  x0 = cut(x, breaks=z$xb, labels=FALSE)
  y0 = cut(y, breaks=z$yb, labels=FALSE)
  return(-as.numeric(z$z[x0, y0]))
}

.trim_traj = function(traj5) {
  traj5$x = round(traj5$x, 2)
  traj5$y = round(traj5$y, 2)
  traj5 = traj5[!duplicated(traj5[, c("x", "y")]), ]
  return(traj5)
}

.getTraj = function(start, N) {
  set.seed(931204)
  out = vector("list", N)
  i = 1
  while(i <= N) {
    gts:::DateStamp(i)
    opt = calibrate(start, fn = find_top, z=z,
                    lower=rep(0.01,2), upper=rep(0.99,2),
                    control=list(REPORT=1, trace=3))
    if(opt$value > -max(z$z)) next
    traj = opt$phases[[1]]$trace$par
    traj = as.data.frame(traj) 
    traj$gen = seq_len(nrow(traj))
    names(traj) = c("x", "y", "gen")
    traj = rbind(c(start, 0), traj)
    out[[i]] = traj
    i = i + 1
  }
  
  xout = do.call(rbind, out)
  # traj0 = aggregate(cbind(x,y) ~ gen, dat=xout, FUN=quantile, prob=0.05)
  # traj0$var = "0"
  traj5 = aggregate(cbind(x,y) ~ gen, dat=xout, FUN=median)
  # traj5$var = "5"
  # traj9 = aggregate(cbind(x,y) ~ gen, dat=xout, FUN=quantile, prob=0.95)
  # traj9$var = "9"
  
  # traj = cbind(traj0, traj5, traj9)
  traj5 = .trim_traj(traj5)
  
  return(traj5)
  
}

.cut_path = function(traj5) {
  if(nrow(traj5)<10) return(traj5)
  x = traj5
  x$x = x$x - tail(x$x,1)
  x$y = x$y - tail(x$y,1)
  x$dist = sqrt(x$x^2 + x$y^2)
  fun = splinefun(x$dist, y=x$gen)
  ndist = exp(seq(from=log(x$dist[6]), to=-3.7, length.out=5))
  ndist = head(ndist, -1)
  ngen = fun(ndist)
  funx = splinefun(x=traj5$gen, y=traj5$x)
  funy = splinefun(x=traj5$gen, y=traj5$y)
  
  xx = c(head(traj5$x, 5), funx(ngen), tail(traj5$x,1))
  yy = c(head(traj5$y, 5), funy(ngen), tail(traj5$y,1))
  out = data.frame(x=xx, y=yy)
  return(out)
}


N = 100
traj5a = .getTraj(start=c(0.01, 0.35), N)
traj5b = .getTraj(start=c(0.01, 0.75), N)
traj5c = .getTraj(start=c(0.25, 0.9), N)
traj5d = .getTraj(start=c(0.4, 0.9), N)
traj5e = .getTraj(start=c(0.5, 0.01), N)
traj5f = .getTraj(start=c(0.6, 0.9), N)
traj5g = .getTraj(start=c(0.9, 0.9), N)
traj5h = .getTraj(start=c(0.99, 0.35), N)

add_line = function(p, data, lwd=1.2) {
  data = .cut_path(data)
  p = p + geom_line(aes(x=x, y=y, z=1), data=data,
                    col="white", alpha=0.5, lwd=lwd)
  p = p + geom_line(aes(x=x, y=y, z=1), data=data,
                    col="red")
  p = p + geom_point(aes(x=x, y=y, z=1), data=data,
                     col="black", cex=0.5)
  return(p)
}

# p <- ggplot(aes(x = x, y = y), data = mt) + geom_point(col="red")
# p <- p + theme_void() + theme_transparent()
 
image.plot(z)
lines(traj5a$x, traj5a$y, col="white", lwd=2)
lines(traj5c$x, traj5c$y, col="white", lwd=2)


p = ggplot(mt, aes(x=x, y=y, z=z))
p = p + geom_contour_filled(show.legend = FALSE)

p = p + theme_void() + theme_transparent()
lwd = 1.2
p = add_line(p, traj5a, lwd=1.5)
# p = add_line(p, traj5b)
p = add_line(p, traj5c, lwd=1.5)
# p = add_line(p, traj5d)
# p = add_line(p, traj5e)
# p = add_line(p, traj5f)
p = add_line(p, traj5g, lwd=1.5)
p = add_line(p, traj5h, lwd=1.5)

# plot(p)

s = sticker(
  # image
  p,
  s_x=1, # slightly to right to appear centered
  s_y=1,
  s_width=2.24,
  s_height=2.24,
  
  # package name
  package="calibrar",
  p_size=25,
  p_color = "white", # 00030A 010101
  p_y = 1.12,
  p_x = 1,
  
  # Output file
  filename="man/figures/logo.png",
  
  # Background colour
  h_fill = "green", # #F0F0F0
  # Border
  # Grey colours: https://www.w3schools.com/colors/colors_shades.asp
  h_color = "#30287F",   # 3F4243 7F2B94 3B2691 4238AF
  h_size = 0,
  # url
  url = "",
  white_around_sticker = TRUE,
  dpi = 300 
)

plot(s)

s = sticker(
  # image
  p,
  s_x=1, # slightly to right to appear centered
  s_y=1,
  s_width=2.24,
  s_height=2.24,
  
  # package name
  package="calibrar",
  p_size=12.5,
  p_color = "white", # 00030A 010101
  p_y = 1.12,
  p_x = 1,
  
  # Output file
  filename="man/figures/logo_small.png",
  
  # Background colour
  h_fill = "green", # #F0F0F0
  # Border
  # Grey colours: https://www.w3schools.com/colors/colors_shades.asp
  h_color = "#30287F",   # 3F4243 7F2B94 3B2691 4238AF
  h_size = 0,
  # url
  url = "",
  white_around_sticker = TRUE,
  dpi = 150 
)
