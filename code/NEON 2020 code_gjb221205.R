library(ggplot2)

# create a vector of site IDs so we can cycle through each site
site = unique(d$Site_ID)

# loop through each site in the vector
for(i in site){
  # create a subset of the data for the current site
  d_sub = d[d$Site_ID==i,]
  
  cd = unique(d_sub$Collection_Date)
  
  for(j in cd){
    #create date subset
    d_subsub = d_sub[d_sub$Collection_Date == j,]
    
    # create a subset of the site data for soils and for plants
    d_soil = d_subsub[d_subsub$Type == "Soil",]  
    d_stem = d_subsub[d_subsub$Type == "Stem",]
    
    # make a vector containing all the soil depth values for the site
    depths = sort(unique(d_soil$Depth_meters))
    
    # create a set of brown colors to match the number of soil depths
    # this uses the hsv color model which you can read about online!
    col.soil = hsv(0.1, seq(0.45, 1, len = length(depths)), 
                   seq(1, 0.4, len = length(depths)))
    
    # print the plots to a file
    png(paste0(i, "_", substr(j, 1, 10), ".png"), width = 6, height = 6, units = "in", res = 600)
    
    # set up the plot area by plotting all data from the site
    plot(d2H ~ d18O, data = d_sub, main = paste(i, j), type = "n")
    
    # now add the soils data, using the colors we created to fill in the symbols
    points(d2H ~ d18O, data = d_soil, pch = 21,
           bg = col.soil[match(d_soil$Depth_meters, depths)])
    
    # now do the same business for the plant data
    
    # make a vector containing all the species names for the site
    # if there are other comments in this data field we may need some clean-up...
    species = (unique(d_stem$Sample_Comments))
    
    # create a set of green colors to match the number of soil depths
    col.stem = hsv(seq(0.25, 0.38, len = length(species)), 
                   seq(0.45, 1, len = length(species)), 
                   seq(1, 0.4, len = length(species)))
    
    # add the stem data
    points(d2H ~ d18O, data = d_stem, pch = 22,
           bg = col.stem[match(d_stem$Sample_Comments, species)])
    
    # add GMWL
    abline(10,8)
    
    # add legend, this might plot over the data in some cases!
    legend("bottomright", c(depths, species), pt.bg = c(col.soil, col.stem), 
           pch = c(rep(21, length(depths)), rep(22, length(species))), 
           ncol = 2, text.width = NA, bty = "n")
    dev.off()    
  }
  

}
