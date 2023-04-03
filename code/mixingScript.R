library(isoWater)
library(assignR)
library(sp)
library(raster)

#Get soil and plant data from wiDB
d = wiDB_data(projects = '00384')
d = d$data

#Subset data to get one site
unique(d$Site_ID)
d.wood = d[d$Site_ID == "WOOD_tower",]

#Get a list of all the bouts by date
cd = unique(d.wood$Collection_Date)

#Precipitation source
isoscape = getIsoscapes("GlobalPrecipMA")
wood.site = SpatialPoints(data.frame(d.wood$Longitude[1], d.wood$Latitude[1]), 
                          crs(isoscape))
##Reality check
plot(isoscape[[1]])
points(wood.site)
##Extract precip values at the site
wood.pcp = extract(isoscape, wood.site)

#Create a list to store all the great results
wood.post = list()

#Loop through each bout and do analysis
for(j in 1:length(cd)){
  #Subset data to get on bout
  d.wood.1 = d.wood[d.wood$Collection_Date == cd[j],]

  #Separate plants and soils
  plants = d.wood.1[d.wood.1$Type == "Stem",]
  soils = d.wood.1[d.wood.1$Type == "Soil",]
  
  #Only run the analysis if we have both plant and soil data!
  if(nrow(plants) * nrow(soils)){
    #Aggregate soil data per depth
    depths = unique(d.wood.1$Depth_meters)
    s1 = soils[soils$Depth_meters >= 0 & soils$Depth_meters < 0.10,]
    s2 = soils[soils$Depth_meters >= 0.1 & soils$Depth_meters < 0.20,]
    s3 = soils[soils$Depth_meters >= 0.2 & soils$Depth_meters < 0.30,]
    s4 = soils[soils$Depth_meters >= 0.3 & soils$Depth_meters < 0.40,]
    s5 = soils[soils$Depth_meters >= 0.4,]
    
    ##Combine all depths into a list for the looping
    s = list(s1, s2, s3, s4, s5)
    
    ##Space to store stats from loop
    sstats = data.frame(d2H = numeric(5), d18O = numeric(5), d2Hsd = numeric(5),
                        d18Osd = numeric(5), HOcov = numeric(5))
    
    ##Loop through each depth and get stats
    for(i in 1:5){
      sstats[i,] = c(mean(s[[i]]$d2H), mean(s[[i]]$d18O), sd(s[[i]]$d2H), 
                     sd(s[[i]]$d18O), cov(s[[i]]$d2H, s[[i]]$d18O))
    }
    
    #Combine Soil and Precip sources
    sstats = rbind(sstats, c(wood.pcp[c(1, 3, 2, 4)], 
                             0.8 * wood.pcp[2] * wood.pcp[4]))
    
    #Make sources into an iso obj
    sources = iso(sstats$d2H, sstats$d18O, sstats$d2Hsd, sstats$d18Osd,
                  sstats$HOcov)
    
    #Define EL slope prior...we should pull from the map on waterisotopes.org
    el = c(6.5, 0.5)
    
    #Create our obs object
    obs = iso(plants$d2H, plants$d18O, rep(2, nrow(plants)), 
              rep(0.5, nrow(plants)), rep(0, nrow(plants)))
    
    #Try the mixing analysis???
    wood.post[[j]] = list()
    for(i in 1:nrow(obs)){
      wood.post[[j]][[i]] = mixSource(obs[i,], sources, el, ngens = 50000, ncores = 3)
    }
  }  

}

#Summary statistics for the first sample
View(wood.post[[1]]$summary)

#Plot of all samples showing surface soil contributions
taxa = unique(plants$Sample_Comments)
par(mar = c(5, 5, 1, 1))
plot(density(wood.post[[1]]$results$s1_fraction), xlim = c(0, 1), ylim = c(0, 7), main = "")
for(i in 2:length(wood.post)){
  lines(density(wood.post[[i]]$results$s1_fraction), col = match(plants$Sample_Comments[i], taxa))
}
legend("topright", legend = taxa, col = seq_along(taxa), lwd = 1)

#Plot of all samples showing deep soil contributions
taxa = unique(plants$Sample_Comments)
par(mar = c(5, 5, 1, 1))
plot(density(wood.post[[1]]$results$s7_fraction), xlim = c(0, 1), ylim = c(0, 7), main = "")
for(i in 2:length(wood.post)){
  lines(density(wood.post[[i]]$results$s7_fraction), col = match(plants$Sample_Comments[i], taxa))
}
legend("topright", legend = taxa, col = seq_along(taxa), lwd = 1)
