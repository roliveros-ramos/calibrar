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


.shadow_text = function(p, x, y, label, family="mono", fontface=1, offset=0.05, size=23, col="black", petit=FALSE) {

  if(petit) size=size/2  
  cols = if(col=="white") c(col, "black") else c(col, "white")
  
  p = p + geom_text(x=x+offset, y=y-offset, label=label, 
                    col=cols[2], size=size, family=family, fontface=fontface)
  p = p + geom_text(x=x, y=y, label=label, col=cols[1], 
                    size=size, family=family, fontface=fontface)
  return(p)
}

add_line = function(p, data, lwd=1.2, opt="blue") {
  # data = .cut_path(data)
  tt2 = data.frame(x = head(data$x,-1), xend = tail(data$x,-1), 
                   y = head(data$y,-1), yend = tail(data$y,-1))
  p = p + geom_path(aes(x=x, y=y, z=1), data=data,
                    col="white", alpha=0.5, lwd=lwd)
  # p = p + geom_line(aes(x=x, y=y, z=1), data=data,
                    # col="red")
  p = p + geom_arrowsegment(aes(x=x, xend=xend, y=y,
                                yend=yend, z=1), data=tt2,
                            col="red", lwd=0.4,
                            arrow_positions = 0.5,
                            arrow_fills = "red",
                            arrows = arrow(type="closed",
                                           length=unit(0.07, units="cm")))
  
  p = p + geom_point(aes(x=x, y=y, z=1), data=head(data, -1),
                     col="black", cex=0.35)
  
  p = p + geom_point(aes(x=x, y=y, z=1), data=tail(data, 1),
                     col="black", fill=opt, cex=0.9, pch=23)

  return(p)
}



make_sticker = function(save = TRUE, behind = TRUE, border = TRUE,
                        u_color = "black", bins = 20, palette = NULL,
                        rev=FALSE, col="black", opt="blue", petit=FALSE) {
  
 
  lab = sprintf("%d: %s", i+1, ifelse(is.null(palette), NA, palette)) 
  
  p = ggplot(mt, aes(x=x, y=y, z=z))
  p = p + geom_contour_filled(show.legend = FALSE, bins=bins) 
  # p = p + scale_fill_brewer(palette = "Spectral")
  if(!is.null(palette)) {
    values = sequential_hcl(bins-1, palette=palette, rev=rev)
    p = p + scale_discrete_manual("fill", values = values)
  }
  
  if(behind) {
    p = .shadow_text(p, x=0.5, y=0.5, label="calibrar", 
                     fontface=2, offset=0.003, size=24, col=col, petit=petit)
  }
  
  # p = add_line(p, traj5a, lwd=1.5)
  # p = add_line(p, traj5b, lwd=1.5)
  p = add_line(p, traj5c, lwd=1.5)
  # p = add_line(p, traj5d, lwd=1.5)
  # p = add_line(p, traj5e, lwd=1.5)
  # p = add_line(p, traj5f, lwd=1.5)
  # p = add_line(p, traj5g, lwd=1.5)
  # p = add_line(p, traj5h, lwd=1.5, opt=opt)
  
  if(!behind) {
    p = .shadow_text(p, x=0.5, y=0.5, label="calibrar", 
                     fontface=2, offset=0.003, size=24, col=col, petit=petit)
  }
  
  p = p + theme_void() + theme_transparent()
  
  # plot(p)
  
  logo = if(!petit) "man/figures/logo.png" else "man/figures/logo_small.png"
    
  s = sticker(
    # image
    p,
    s_x=1, # slightly to right to appear centered
    s_y=1,
    s_width=2.24,
    s_height=2.24,
    # package name
    package="",
    p_size=25,
    p_color = "white", # 00030A 010101
    p_y = 1.12,
    p_x = 1,
    # Output file
    filename=ifelse(save, sprintf("man/figures/%03dlogo.png", i <<- i+1), logo),
    # Background colour
    h_fill = "green", # #F0F0F0
    # Border
    # Grey colours: https://www.w3schools.com/colors/colors_shades.asp
    h_color = "black",   # 3F4243 7F2B94 3B2691 4238AF
    h_size = ifelse(border, 2, 0),
    # url
    url = ifelse(petit, 
                 "Parameter estimation for complex models",
                 "  Parameter estimation for complex models"),
    u_size = ifelse(petit, 1.91, 3.82),
    u_color = u_color,
    white_around_sticker = TRUE,
    dpi = ifelse(petit, 150, 300) 
  )
  
  if(save) plot(s)
  mtext(lab, 3, adj=0.05)
  
  return(invisible(s))
  
}




