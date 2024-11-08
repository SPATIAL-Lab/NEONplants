violin.prep = function(ss){
  ssdf = data.frame("Species" = character(), "Depth" = numeric())
  for(i in seq_along(ss$taxa)){
    ssdf = rbind(ssdf, data.frame("Species" = rep(ss$taxa[i]),
                                  "Depth" = ss$source.depths[i,]))
  }
  return(ssdf)
}

violin = function(d, ...){
  par(mar = c(8, 5, 3, 1))
  cats = unique(d[, 1])
  ylim = c(5, 1)
  xlim = c(0.5, length(cats) + 0.5)
  
  plot(0, 0, xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
       xlab = "", ylab = "", ...)
  
  for(i in seq_along(cats)){
    d.sub = d[d$Species == cats[i], 2]
    d.dens = density(d.sub)
    xscale = 0.45 / max(d.dens$y)
    
    polygon(c(d.dens$y * xscale + i, i - rev(d.dens$y) * xscale), 
            c(d.dens$x, rev(d.dens$x)), ...)
  }
  
  axis(1, at = seq(length(cats)), labels = FALSE, lwd = 0, lwd.ticks = 1)
  axis(2, at = seq(1:5), labels = c("0-10 cm", "10-20 cm", "20-40 cm",
                                    ">40 cm", "GW"), las = 1)
  
  for(i in seq_along(cats)){
    text(i, 5.4, bquote(italic(.(cats[i]))), 
         srt = 45, xpd = NA, adj = 1)
  }
}
