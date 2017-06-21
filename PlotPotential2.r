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


PlotPotential_list = function (res_list, title = "", xlab.text, ylab.text, cutoff = 0.5, 
                           plot.contours = TRUE, binwidth = 0.2, bins = NULL) 
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
  
  lapply(res_list, res_to_df)
  
  
  
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

