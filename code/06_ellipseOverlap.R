library(SIBER)
source("code/multiLikeOverlap.R")

# Using 'raw' data ----
## Read data
p = read.csv("data/plants.csv")

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
  
  ## Limit to at least 4 species
  if(length(spi) > 3){
    pss = pss[pss$group %in% spi,]
    
    ## SIBER object
    ps.t = createSiberObject(pss)
    
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
    mlo[[i]] = NULL
  }
}

## Names and counts
names(mlo) = bouts
spnum = unlist(lapply(mlo, length)) - 2
spnum = spnum[spnum > 0]

## Pull out overlap values and combine
overlaps = unlist(lapply(mlo, tail, n = 1))
overlaps = data.frame("Bout" = names(spnum), "SpeciesNum" = spnum, 
                      "Overlap" = overlaps)

## Plot and write
plot(density(overlaps$Overlap))
write.csv(overlaps, "out/MLO_raw.csv")

# Using source water data ----
