library(SIBER)
source("code/multiLikeOverlap.R")

# Using 'raw' data ----
## Read data
p = read.csv("out/plantsNoLow.csv")

## Get bits for SIBER
ps = data.frame("iso1" = p$d2H, "iso2" = p$d18O, "group" = p$Species, 
                "community" = p$Bout)

## Make sequential values for community
bouts = unique(ps$community)
ps$community = match(ps$community, bouts)

## Space for output
mlo = list()

## Plotting args
group.ellipses.args = list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)

## Loop through bouts
for(i in seq_along(bouts)){
  ## Pull out bout
  pss = ps[ps$community == i,]
  spf = table(pss$group)
  spi = names(spf)[spf > 3]
  
  pss = pss[pss$group %in% spi,]

  if(length(unique(pss$group)) > 1){
    ## SIBER object
    ps.t = suppressWarnings(createSiberObject(pss))
    
    ## Plot for fun
    par(mar = c(5,5,3,1))
    plotSiberObject(ps.t,                   
                    ax.pad = 5, 
                    hulls = FALSE, community.hulls.args, 
                    ellipses = TRUE, group.ellipses.args,
                    group.hulls = FALSE, group.hull.args,
                    bty = "L",
                    iso.order = c(2, 1),
                    xlab = expression(delta^18*O~'\u2030'),
                    ylab = expression(delta^2*H~'\u2030'))
    box()
    legend("topleft", ps.t$all.groups, col = seq(length(ps.t$all.groups)), lty = 1,
           bty = "n")
    
    ## Run MLE using my multi version...so fast
    mlo[[i]] = multiLikOverlap(ps.t$all.groups, ps.t, n = 500)
    
    title(paste(bouts[i], ":", round(mlo[[i]]["Fraction Overlap"], 2)))
  }else{
    mlo[[i]] = NA
  }
}

## Names and counts
names(mlo) = bouts
spnum = unlist(lapply(mlo, length)) - 1

## Pull out overlap values and combine
overlaps = unlist(lapply(mlo, tail, n = 1))
overlaps = data.frame("Bout" = names(spnum), "SpeciesNum" = spnum, 
                      "Overlap.raw" = overlaps)

# Using source water data, all posterior draws ----
## Get bout names and matching sites
p = read.csv("out/plantsNoLow.csv")
bouts = data.frame("Bout" = unique(p$Bout))
bouts$Site = p$Site_ID[match(bouts$Bout, p$Bout)]

## Space for output
mlo = list()

## Plotting args
group.ellipses.args = list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)

for(i in seq_along(bouts$Bout)){
  ## Load and prep data for bout
  load(paste0("out/", bouts$Site[i], "/", bouts$Bout[i], ".rda"))
  ps = data.frame("iso1" = numeric(), "iso2" = numeric(), "group" = character(), 
                  "community" = numeric())
  for(j in seq_along(smix)){
    psa = data.frame("iso1" = smix[[j]]$results$mixture_d2H, 
                     "iso2" = smix[[j]]$results$mixture_d18O, 
                     "group" = rep(names(smix)[j]), 
                     "community" = rep(1))
    ps = rbind(ps, psa)
  }
  
  if(length(unique(ps$group)) > 1){
    ## SIBER object
    ps.t = createSiberObject(ps)
    
    ## Plot for fun
    par(mar = c(5,5,3,1))
    plotSiberObject(ps.t,                   
                    ax.pad = 5, 
                    hulls = FALSE, community.hulls.args, 
                    ellipses = TRUE, group.ellipses.args,
                    group.hulls = FALSE, group.hull.args,
                    bty = "L",
                    iso.order = c(2, 1),
                    xlab = expression(delta^18*O~'\u2030'),
                    ylab = expression(delta^2*H~'\u2030'))
    box()
    legend("topleft", ps.t$all.groups, col = seq(length(ps.t$all.groups)), lty = 1,
           bty = "n")
    
    ## Run MLE using my multi version...so fast
    mlo[[i]] = multiLikOverlap(ps.t$all.groups, ps.t, n = 500)
    
    title(paste(bouts$Bout[i], ":", round(mlo[[i]]["Fraction Overlap"], 2)))
  }else{
    mlo[[i]] = NA
  }
}

overlap.add = unlist(lapply(mlo, tail, n = 1))
overlaps = cbind(overlaps, "Overlap.post" = overlap.add)

# Using source water data, medians ----
for(i in seq_along(bouts$Bout)){
  ## Load and prep data for bout
  load(paste0("out/", bouts$Site[i], "/", bouts$Bout[i], ".rda"))
  ps = data.frame("iso1" = numeric(), "iso2" = numeric(), "group" = character(), 
                  "community" = numeric())
  for(j in seq_along(smix)){
    psa = data.frame("iso1" = smix[[j]]$summary["mixture_d2H", "50%"], 
                     "iso2" = smix[[j]]$summary["mixture_d18O", "50%"], 
                     "group" = names(smix)[j], 
                     "community" = 1)
    ps = rbind(ps, psa)
  }
  
  spf = table(ps$group)
  spi = names(spf)[spf > 3]
  
  ps = ps[ps$group %in% spi,]
  
  if(length(unique(ps$group)) > 1){
    ## SIBER object
    ps.t = suppressWarnings(createSiberObject(ps))
    
    ## Plot for fun
    par(mar = c(5,5,3,1))
    plotSiberObject(ps.t,                   
                    ax.pad = 5, 
                    hulls = FALSE, community.hulls.args, 
                    ellipses = TRUE, group.ellipses.args,
                    group.hulls = FALSE, group.hull.args,
                    bty = "L",
                    iso.order = c(2, 1),
                    xlab = expression(delta^18*O~'\u2030'),
                    ylab = expression(delta^2*H~'\u2030'))
    box()
    legend("topleft", ps.t$all.groups, col = seq(length(ps.t$all.groups)), lty = 1,
           bty = "n")
    
    ## Run MLE using my multi version...so fast
    mlo[[i]] = multiLikOverlap(ps.t$all.groups, ps.t, n = 500)
    
    title(paste(bouts$Bout[i], ":", round(mlo[[i]]["Fraction Overlap"], 2)))
  }else{
    mlo[[i]] = NA
  }
}

overlap.add = unlist(lapply(mlo, tail, n = 1))
overlaps = cbind(overlaps, "Overlap.postMed" = overlap.add)

## Plot and write
plot(density(overlaps$Overlap.raw, na.rm = TRUE))
lines(density(overlaps$Overlap.postMed, na.rm = TRUE), col = "red")
lines(density(overlaps$Overlap.post, na.rm = TRUE), col = "blue")

plot(overlaps$Overlap.raw, overlaps$Overlap.post, pch = 21, 
     bg = overlaps$SpeciesNum)
plot(overlaps$Overlap.raw, overlaps$Overlap.postMed, pch = 21, 
     bg = overlaps$SpeciesNum)

write.csv(overlaps, "out/MLO.csv", row.names = FALSE)
