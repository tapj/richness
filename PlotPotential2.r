PlotPotential2 = function (res, title = "", xlab.text, ylab.text, cutoff = 0.5, 
          plot.contours = TRUE, binwidth = 0.2, bins = NULL) 
{
  cut.potential <- max(apply(res$pots, 1, min)) + cutoff * 
    abs(max(apply(res$pots, 1, min)))
  pots <- res$pots
  pots[pots > cut.potential] <- cut.potential
  intp <- tgp::interp.loess(as.vector(res$pars), as.vector(res$xis), 
                            as.vector(pots), gridlen = 2 * dim(pots))
  xy <- expand.grid(intp$x, intp$y)
  z <- as.vector(intp$z)
  z[is.na(z)] <- max(na.omit(z))
  potential <- NULL
  df <- data.frame(list(bg.var = xy[, 1], phylotype = xy[, 
                                                         2], potential = z))
  bg.var <- NULL
  phylotype <- NULL
  p <- ggplot2::ggplot(df, aes(bg.var, phylotype, z = potential)) + 
    geom_tile(aes(fill = potential)) + scale_fill_gradientn(colours = topo.colors(10))
  if (plot.contours) {
    if (!is.null(bins)) {
      warning("bins argument is overriding the binwidth argument!")
      p <- p + stat_contour(bins = bins)
    }
    else {
      p <- p + stat_contour(binwidth = binwidth)
    }
  }
  p <- p + xlab(xlab.text) + ylab(ylab.text) + labs(title = title)
  p
}


Potential_list = function (res_list, title = "", xlab.text, ylab.text, cutoff = 0.5, 
                           plot.contours = FALSE, binwidth = 0.2, bins = NULL, output="ani.gif") 
{
  
  res_to_df = function(res){
  cut.potential <- max(apply(res$pots, 1, min)) + cutoff * 
    abs(max(apply(res$pots, 1, min)))
  pots <- res$pots
  pots[pots > cut.potential] <- cut.potential
  intp <- tgp::interp.loess(as.vector(res$pars), as.vector(res$xis), 
                            as.vector(pots), gridlen = 2 * dim(pots))
  xy <- expand.grid(intp$x, intp$y)
  z <- as.vector(intp$z)
  z[is.na(z)] <- max(na.omit(z))
  potential <- NULL
  df <- data.frame(list(bg.var = xy[, 1], phylotype = xy[, 
                                                         2], potential = z))
  bg.var <- NULL
  phylotype <- NULL
  
  return(df)
  }
  
  df_list = lapply(res_list, res_to_df)
  
  tf = tweenr::tween_states(df_list, 
                            tweenlength=2, statelength=1, 
                            ease="linear", nframes = 300)
  
  return(tf)
  
}  


Potential_animate=function(tf, output = "microbiota_animation.gif", speed = 10) {

animation::ani.options(interval = 1/speed)

animation::saveGIF({   
 for(i in 1:length(table(tf[,4]))) {
 p <- ggplot2::ggplot(tf[tf[,".frame"] %in% i,], aes(bg.var, phylotype, z = potential, frame = .frame)) 
  p <- p + geom_tile(aes(fill = potential)) + scale_fill_gradientn(colours = topo.colors(10))
  if (plot.contours) {
    if (!is.null(bins)) {
      warning("bins argument is overriding the binwidth argument!")
      p <- p + stat_contour(bins = bins)
    }
    else {
      p <- p + stat_contour(binwidth = binwidth)
    }
  }
  p <- p + xlab(xlab.text) + ylab(ylab.text) + labs(title = title) + ylim(-0.2, 0.2) + guides(fill=FALSE)
  print(p)
  #output = paste0("plot",i,".png")
  #ggsave(plot = p, filename = output)
  
 }
}, movie.name = output, clean = TRUE)
  
  
# gganimate did not work so I had to use saveGIF instead
  #gganimate(p, output, title_frame = F, ani.width = 400, 
  #          ani.height = 400)
  
  #gganimate(p, output, title_frame = F) 
  
}

