# Prep ----
## Read data
p = read.csv("out/plantsNoLow.csv")
s = read.csv("data/soils.csv")

## Bouts
bouts = sort(unique(c(p$Bout, s$Bout)))

# Plot ----
for(i in bouts){
  ## Bout subset
  p_sub = p[p$Bout == i,]
  s_sub = s[s$Bout == i,]
  xlim = range(c(p_sub$d18O.oc, s_sub$d18O))
  ylim = range(c(p_sub$d2H, s_sub$d2H))
  
  ## Depths and species
  species = (unique(p_sub$Species))
  depths = sort(unique(s_sub$Depth))
  
  ## Create a set of brown colors to match the number of soil depths
  col.soil = hsv(0.1, seq(0.45, 1, len = length(depths)), 
                 seq(1, 0.4, len = length(depths)))
  
  ## create a set of green colors to match the number of species
  col.stem = hsv(seq(0.25, 0.38, len = length(species)), 
                 seq(0.45, 1, len = length(species)), 
                 seq(1, 0.4, len = length(species)))
  
  ## Output device
  png(paste0("out/", i, ".png"), width = 6, height = 6, units = "in", res = 600)
  par(mar = c(5, 5, 3, 3))
  
  ## Set up the plot area
  plot(d2H ~ d18O.oc, data = p_sub, xlim = xlim, ylim = ylim,
       main = paste(i, ":", format(as.Date(p_sub$Collection_Date[1]), "%Y-%m")), 
       type = "n", xlab = expression(delta^{18}*"O (VSMOW)"),
       ylab = expression(delta^2*"H (VSMOW)"))
  
  ## Add the soils data
  points(d2H ~ d18O, data = s_sub, pch = 21,
         bg = col.soil[match(s_sub$Depth, depths)])
  
  ## Add the stem data
  points(d2H ~ d18O.oc, data = p_sub, pch = 22,
         bg = col.stem[match(p_sub$Species, species)])
  
  ## Add GMWL
  abline(10,8)
  
  ## Add legend, this might plot over the data in some cases!
  legend("bottomright", c(species, depths), pt.bg = c(col.stem, col.soil), 
         pch = c(rep(22, length(species)), rep(21, length(depths))), 
         ncol = 2, text.width = NA, bty = "n")
  dev.off()    
}
