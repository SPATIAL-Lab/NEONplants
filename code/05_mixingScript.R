library(isoWater)
library(assignR)
library(terra)
library(ggplot2)

# Prep ----
## Read data
p = read.csv("out/plantsNoLow.csv")
s = read.csv("data/soils.csv")

## GW source map
isoscape = getIsoscapes("USGround")

## Sites
sites = unique(data.frame("ID" = p$Site_ID, "lon" = p$Longitude, 
                          "lat" = p$Latitude))
sites = vect(sites, crs = "WGS84")
sites = project(sites, isoscape)

## GW at sites
gw = extract(isoscape, sites, method = "bilinear")

## Space for summary stats
mixes = list()

## Space to track whether plant values are bounded by sources
pbound = p[c("Site_ID", "Sample_ID", "Bout", "Species")]
pbound$d2H.bound = pbound$d18O.bound = rep(NA)

# Mix source ----
## Loop through each site and bout and do analysis
for(i in 1:length(sites)){
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
    snames = c(depths, "Groundwater")

    ## Make sources into an iso obj
    sources = iso(data.frame(svals))

    ## EL slope and evap prior parameters
    el = c(2.5, 0.5)
    eprior = c(0.2, 1)
    
    ## Iso object for plant samples
    obs = iso(p.bout$d2H, p.bout$d18O.oc, rep(2, nrow(p.bout)), 
              rep(0.5, nrow(p.bout)), rep(0, nrow(p.bout)))
    
    ## Mixing analysis
    smix = list()
    for(k in 1:nrow(obs)){
      smix[[k]] = mixSource(obs[k,], sources, el, edist = "gamma",
                            eprior = eprior, ngens = 5e4, ncores = 3)
      names(smix)[k] = p.bout$Species[k]
      names(smix[[k]]$results)[3:(2 + length(snames))] = snames
      row.names(smix[[k]]$summary)[4:(3 + length(snames))] = snames
      
      ## Compile summaries
      mix = list("Bout" = bid, "Species" = p.bout$Species[k], 
                 "Mix" = smix[[k]]$summary)
      mixes[[length(mixes) + 1]] = mix
      
      ## Populate pbound
      if(obs$H[k] >= min(sources$H - 2 * sources$Hsd) & 
         obs$H[k] <= max(sources$H + 2 * sources$Hsd)){
        pbound$d2H.bound[pbound$Sample_ID == p.bout$Sample_ID[k]] = TRUE
      } else{
        pbound$d2H.bound[pbound$Sample_ID == p.bout$Sample_ID[k]] = FALSE
      }
      
      if(obs$O[k] >= min(sources$O - 2 * sources$Osd) & 
         obs$O[k] <= max(sources$O + 2 * sources$Osd)){
        pbound$d18O.bound[pbound$Sample_ID == p.bout$Sample_ID[k]] = TRUE
      } else{
        pbound$d18O.bound[pbound$Sample_ID == p.bout$Sample_ID[k]] = FALSE
      }
    }
    
    ## Write results
    sources = cbind("Depth" = snames, sources)
    save(smix, sources, file = file.path("out", sid, paste0(bid, ".rda")))
    
  }
}

## Summarize pbound and save
sum(pbound$d18O.bound) / nrow(pbound)
sum(pbound$d2H.bound) / nrow(pbound)

# Parse mixing summaries ----
nr = length(mixes)
fn = row.names(mixes[[1]]$Mix)[-c(2, 3,9, 10)]
mix = data.frame("Bout" = character(nr), "Species" = character(nr), 
                 "E.med" = numeric(nr), "E.IQR" = numeric(nr),
                 "X0.10cm.med" = numeric(nr), "X0.10cm.IQR" = numeric(nr),
                 "X10.20cm.med" = numeric(nr), "X10.20cm.IQR" = numeric(nr),
                 "X20.40cm.med" = numeric(nr), "X20.40cm.IQR" = numeric(nr),
                 "X40.cm.med" = numeric(nr), "X40.cm.IQR" = numeric(nr),
                 "GW.med" = numeric(nr), "GW.IQR" = numeric(nr))

for(i in seq_along(mixes)){
  mix$Bout[i] = mixes[[i]]$Bout
  mix$Species[i] = mixes[[i]]$Species
  for(j in seq_along(fn)){
    k = match(fn[j], row.names(mixes[[i]]$Mix))
    if(is.na(k)){
      mix[i, 3 + (j - 1) * 2] = mix[i, 4 + (j - 1) * 2] = NA
    }else{
      mix[i, 3 + (j - 1) * 2] = mixes[[i]]$Mix[k, 5]
      mix[i, 4 + (j - 1) * 2] = mixes[[i]]$Mix[k, 6] - mixes[[i]]$Mix[k, 4]
    }
  }
}

write.csv(mix, "out/mixStats.csv", row.names = FALSE)

# Plot results: density by species ----
## Function
plot.post = function(smix, bid){
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
    plot(d[[1]], main = paste(bid, sp.all[i], sep = ":"), xlab = "", xlim = c(0, 1), lwd = 2,
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
          d[[j]] = density(smix[[sp.ind[k]]]$results[, 2 + j], from = 0, to = 1,
                           bw = 0.05)
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
}

## Loop
for(i in 1:length(sites)){
  sid = sites$ID[i]
  bouts = sort(unique(p[p$Site_ID == sid,]$Bout))
  for(j in seq_along(bouts)){
    load(file.path("out", sid, paste0(bouts[j], ".rda")))
    plot.post(smix, bouts[j])
  }
}

# Plot results: violin by bout ----
## Function
plot.sources = function(smix, bid){
  
  r = smix[[1]]$results
  r$species = rep(names(smix)[1])
  for(i in seq_along(smix)[-1]){
    ra = smix[[i]]$results
    ra$species = rep(names(smix)[i])
    r = rbind(r, ra)
  }
  
  rv = data.frame("species" = character(), "source" = character(), "fraction" = numeric())
  for(i in 1:(ncol(r) - 5)){
    rva = data.frame("species" = r$species, "fraction" = r[, 2 + i])
    rva$source = rep(names(r)[i + 2])
    rv = rbind(rv, rva)
  }
  
  png(paste0("out/sources/", bid, ".png"), width = 9, height = 4, units = "in", res = 600)
  print(ggplot(rv, aes(x = species, y = fraction, fill = source)) +
    geom_violin() + 
    labs(title = bid, x = "Species", y = "Fraction", fill = "Sources"))
  dev.off()
  
}

## Loop
for(i in 1:length(sites)){
  sid = sites$ID[i]
  bouts = sort(unique(p[p$Site_ID == sid,]$Bout))
  for(j in seq_along(bouts)){
    load(file.path("out", sid, paste0(bouts[j], ".rda")))
    plot.sources(smix, bouts[j])
  }
}



