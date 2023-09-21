
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
library(ggarchery)
library(colorspace)

source("man/figures/sticker_functions.R")

mt = expand.grid(x=seq_len(ncol(volcano)), 
                 y=seq_len(nrow(volcano)))
mt$z = as.numeric(t(volcano))
mt$x = (mt$x-1.01)/(max(mt$x)-1)
mt$y = (mt$y-1.01)/(max(mt$y)-1)
xout = seq(0.01, 0.99,length.out=1000)

z = interpolate(mt$x, mt$y, mt$z, xout=xout, yout=xout)

mt2 = expand.grid(x=xout, y=xout)
mt2$z = as.numeric(z$z)

z$xb = gts:::.getBreaks(z$x)
z$yb = gts:::.getBreaks(z$y)

N = 100
set.seed(931204)
traj5a = .getTraj(start=c(0.01, 0.35), N)
traj5b = .getTraj(start=c(0.01, 0.75), N)
traj5c = .getTraj(start=c(0.25, 0.9), N)
traj5d = .getTraj(start=c(0.4, 0.9), N)
traj5e = .getTraj(start=c(0.5, 0.01), N)
traj5f = .getTraj(start=c(0.6, 0.99), N)
traj5g = .getTraj(start=c(0.9, 0.9), N)
traj5h = .getTraj(start=c(0.99, 0.35), N)

i = 0

par(mfrow=c(3,4), mar=c(1,1,1,1), oma=c(1,1,1,1))

make_sticker(palette="GnBu")
make_sticker(palette="GnBu", rev=TRUE)
make_sticker(palette="PuBu", rev=FALSE) 
make_sticker(palette="Plasma", rev=FALSE) 
make_sticker(palette="Viridis", rev=FALSE) 
make_sticker(palette="Inferno", rev=FALSE, opt="purple4") 
make_sticker(palette="Inferno", rev=FALSE, col="white", opt="purple4") 
make_sticker(palette="YlGnBu", rev=FALSE) 
make_sticker(palette="BuPu", rev=FALSE) 
make_sticker(palette="ag_Sunset", rev=FALSE) 
make_sticker(palette="Lajolla", rev=FALSE) 
make_sticker(palette="SunsetDark", rev=FALSE) 

make_sticker(palette="Inferno", rev=FALSE, col="white", opt="purple4", save=FALSE) 
make_sticker(palette="Inferno", rev=FALSE, col="white", opt="purple4", save=FALSE, petit=TRUE) 

beepr::beep(4)

xx = .cut_path(traj5c)
