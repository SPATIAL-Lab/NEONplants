library(isoWater)
library(assignR)
library(terra)

# Prep ----
## Read data
p = read.csv("out/plantsNoLow.csv")
s = read.csv("data/soils.csv")

## GW source map
isoscape = getIsoscapes("USGround")

## Sites
sites = unique(data.frame("ID" = p$Site_ID, "lon" = p$Longitude, "lat" = p$Latitude))
sites = vect(sites, crs = "WGS84")
sites = project(sites, isoscape)

## GW at sites
gw = extract(isoscape, sites, method = "bilinear")

# Mix source ----
## Loop through each site and bout and do analysis
for(i in 10:length(sites)){
  ## Subset for site
  sid = sites$ID[i]
  p.site = p[p$Site_ID == sid,]
  s.site = s[s$Site_ID == sid,]
  
  ## Storage space
  if(!(dir.exists(file.path("out", sid)))){
    dir.create(file.path("out", sid))
  }
  
  ## List of bouts
  bouts = unique(p.site$Bout)
  
  for(j in seq_along(bouts)){
    ## Subset for bout
    bid = bouts[j]
    p.bout = p.site[p.site$Bout == bid,]
    s.bout = s.site[s.site$Bout == bid,]
    
    ## Prep soil sources
    depths = unique(s.bout$Depth)
    
    ## Space to store soil stats
    svals = matrix(ncol = 5, nrow = 0)
    
    ## Loop through each depth and get stats
    for(k in seq_along(depths)){
      s.depth = s.bout[s.bout$Depth == depths[k],]
      ## Only calculate covariance if n > 3, otherwise assume R = 0.8
      if(nrow(s.depth) > 2){
        svals = rbind(svals, c(mean(s.depth$d2H), mean(s.depth$d18O), 
                               sd(s.depth$d2H), sd(s.depth$d18O), 
                               cov(s.depth$d2H, s.depth$d18O)))
      } else{
        svals = rbind(svals, c(mean(s.depth$d2H), mean(s.depth$d18O), 
                               sd(s.depth$d2H), sd(s.depth$d18O), 
                               0.8 * sd(s.depth$d2H) * sd(s.depth$d18O)))
      }
    }
    
    ## Fill in missing values
    for(k in 3:4){
      svals[is.na(svals[,k]), k] = mean(svals[,k], na.rm = TRUE)
    }
    for(k in seq_along(svals[,1])){
      if(is.na(svals[k, 5])){
        svals[k, 5] = svals[k, 3] * svals[k, 4] * 0.8
      }
    }

    ## Add groundwater assume R = 0.8
    svals = rbind(svals, c(gw$`d2h_1-10m`[i], gw$`d18o_1-10m`[i],
                           gw$`d2h_sd_1-10m`[i], gw$`d18o_sd_1-10m`[i],
                           0.8 * gw$`d2h_sd_1-10m`[i] * gw$`d18o_sd_1-10m`[i]))

    ## Make sources into an iso obj
    sources = iso(data.frame(svals))

    ## EL slope and evap prior parameters
    el = c(2.5, 0.5)
    eprior = c(0.2, 1)
    
    ## Iso object for plant samples
    obs = iso(p.bout$d2H, p.bout$d18O, rep(2, nrow(p.bout)), 
              rep(0.5, nrow(p.bout)), rep(0, nrow(p.bout)))
    
    ## Mixing analysis
    smix = list()
    for(k in 1:nrow(obs)){
      smix[[k]] = mixSource(obs[k,], sources, el, edist = "gamma",
                            eprior = eprior, ngens = 5e4, ncores = 3)
      names(smix)[k] = p.bout$Species[k]
    }
    
    ## Write results
    save(smix, file = file.path("out", sid, paste0(bid, ".rda")))
  }
}

# Plot results ----
species = names(smix)
sp.all = unique(species)
sp.ind = match(species, sp.all)

## Cycle through species
for(i in seq_along(sp.all)){
  sp.ind = which(species == sp.all[i])
  
  ## Space for density
  d = list()
  d.all = numeric()
  for(j in 1:5){
    d[[j]] = density(smix[[sp.ind[1]]]$results[, 2 + j], from = 0, to = 1)
    d.all = append(d.all, d[[j]]$y)
  }
  
  ## Plot first source
  plot(d[[1]], main = sp.all[i], xlab = "", xlim = c(0, 1), lwd = 2,
       ylim = c(0, max(d.all)))
  
  ## Add others
  for(j in 2:5){
    lines(d[[j]], col = j, lwd = 2)
  }
  
  if(length(sp.ind) > 1){
    for(k in 2:length(sp.ind)){
      d = list()
      d.all = numeric()
      for(j in 1:5){
        d[[j]] = density(smix[[sp.ind[k]]]$results[, 2 + j], from = 0, to = 1)
        d.all = append(d.all, d[[j]]$y)
      }
      
      ## Add sources
      for(j in 1:5){
        lines(d[[j]], col = j, lty = k, lwd = 2)
      }
    }
  }
  
  legend("topright", c("0-10 cm", "10-20 cm", "20-40 cm", "40+ cm", "Groundwater"),
         col = 1:5, lty = 1, lwd = 2)
}



# Code residue, some may be useful ----    
    
    #Append names to the samples
    names(wood.post[[i]]) = plants$Sample_ID
    names(wood.post)[i] = cd[i]

    #Plot data
    taxa = unique(plants$Sample_Comments)
    plot(plants$d18O, plants$d2H, xlab = "d18O", ylab = "d2H",
         xlim = range(c(plants$d18O, soils$d18O)),
         ylim = range(c(plants$d2H, soils$d2H)))
    points(sstats$d18O, sstats$d2H, pch = 21, bg = seq(nrow(sstats)))
    for(j in 2:length(taxa)){
      points(plants$d18O[plants$Sample_Comments == taxa[j]],
             plants$d2H[plants$Sample_Comments == taxa[j]], col = j)
    }
    abline(10, 8)
    legend("bottomright", legend = taxa, col = seq_along(taxa), pch = 1)

    #Plots of each source contribution
    par(mar = c(5, 5, 4, 1))
    for(j in 1:nrow(sources)){
      plot(density(smix[[1]]$results[, 2+j]), 
           xlim = c(0, 1), ylim = c(0, 7), 
           main = paste("source", j))
      for(k in 2:length(smix)){
        lines(density(smix[[k]]$results[, 2+j]), 
              col = match(names(smix)[k], unique(names(smix))))
      }
      legend("topright", legend = unique(names(smix)), 
             col = seq_along(unique(names(smix))), lwd = 1)
    }
    


#Print out effective sample sizes for review
for(i in 1:length(wood.post)){
  for(j in 1:length(wood.post[[i]])){
    print(names(wood.post)[i])
    print(names(wood.post[[i]])[j])
    print(wood.post[[i]][[j]]$summary[,9])
  }
}

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
