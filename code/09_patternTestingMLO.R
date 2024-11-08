library(terra)

## Water source info
mlo = read.csv("out/MLO.csv")

## Properties
bouts = read.csv("out/Bouts.csv")
sites = values(vect("out/Sites.shp"))

## Work with bouts w/ 4 or more taxa
mlo = mlo[mlo$SpeciesNum > 2,]

# Wetter sites have less diverse water use ----
## Add site means to sites
mlo$Site = substr(mlo$Bout, 1, 4)
sites$Overlap.raw = sites$Overlap.post = sites$Overlap.postMed = rep(0)
for(i in seq_along(sites$Site)){
  sites$Overlap.raw[i] = 
    mean(mlo$Overlap.raw[mlo$Site == sites$Site[i]])
  sites$Overlap.post[i] = 
    mean(mlo$Overlap.post[mlo$Site == sites$Site[i]])
  sites$Overlap.postMed[i] = 
    mean(mlo$Overlap.postMed[mlo$Site == sites$Site[i]])
}

## Add all bouts to mlo.sites
mlo.sites = merge(mlo, sites)
boxplot(Overlap.raw ~ Site, mlo)
boxplot(Overlap.post ~ Site, mlo)
boxplot(Overlap.postMed ~ Site, mlo)

## Based on MAP - no
plot(sites$MAP, sites$Overlap.postMed)
summary(lm(Overlap.postMed ~ MAP, sites))

## Based on VPD - yes
plot(sites$VPDmax, sites$Overlap.postMed) ####
summary(lm(Overlap.postMed ~ VPDmax, sites))

# Wetter bouts have less diverse water use ----
## Add bouts to mlo.bouts
mlo.bouts = merge(mlo, bouts, by = "Bout")

## Based on 10 day precip
plot(mlo.bouts$Precip10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$Precip10, mlo.bouts$Overlap.postMed)

## Based on 30 day precip
plot(mlo.bouts$Precip30, mlo.bouts$Overlap.raw)
plot(mlo.bouts$Precip30, mlo.bouts$Overlap.postMed)

## Based on 60 day precip
plot(mlo.bouts$Precip60, mlo.bouts$Overlap.raw)
plot(mlo.bouts$Precip60, mlo.bouts$Overlap.postMed)

## Based on 1 day shallow SWC
plot(mlo.bouts$VWC1.1, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC1.1, mlo.bouts$Overlap.postMed) ####

## Based on 1 day mid SWC
plot(mlo.bouts$VWC3.1, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC3.1, mlo.bouts$Overlap.postMed) ####

## Based on 1 day deep SWC
plot(mlo.bouts$VWC5.1, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC5.1, mlo.bouts$Overlap.postMed) ####

## Based on 10 day shallow SWC
plot(mlo.bouts$VWC1.10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC1.10, mlo.bouts$Overlap.postMed)

## Based on 10 day mid SWC
plot(mlo.bouts$VWC3.10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC3.10, mlo.bouts$Overlap.postMed)

## Based on 10 day deep SWC
plot(mlo.bouts$VWC5.10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC5.10, mlo.bouts$Overlap.postMed)

# Bouts with more diverse water use have higher LH flux ----
plot(mlo.bouts$lh10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$lh10, mlo.bouts$Overlap.postMed)

# Bouts with more diverse water use have lower LH variation ----
plot(mlo.bouts$lhVar, mlo.bouts$Overlap.raw)
plot(mlo.bouts$lhVar, mlo.bouts$Overlap.postMed)

# Wetter bouts have less diverse water use than local ave ----
## Calculate bout anomalies
mlo.bAnom = mlo.bouts[0, ]
for(i in sites$Site){
  mlo.sub = mlo.bouts[mlo.bouts$Site.x == i,]
  for(j in c(3:5, 9:19)){
    mlo.sub[, j] = mlo.sub[, j] - mean(mlo.sub[, j], na.rm = TRUE)
  }
  mlo.bAnom = rbind(mlo.bAnom, mlo.sub)
}

## Based on 10 day precip
plot(mlo.bAnom$Precip10, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$Precip10, mlo.bAnom$Overlap.postMed)

## Based on 30 day precip
plot(mlo.bAnom$Precip30, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$Precip30, mlo.bAnom$Overlap.postMed)

## Based on 60 day precip
plot(mlo.bAnom$Precip60, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$Precip60, mlo.bAnom$Overlap.postMed) ####
summary(lm(Overlap.postMed ~ Precip60, mlo.bAnom))

## Based on 1 day shallow SWC
plot(mlo.bAnom$VWC1.1, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$VWC1.1, mlo.bAnom$Overlap.postMed) ####
summary(lm(Overlap.postMed ~ VWC1.1, mlo.bAnom))

## Based on 1 day mid SWC
plot(mlo.bAnom$VWC3.1, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$VWC3.1, mlo.bAnom$Overlap.postMed)

## Based on 1 day deep SWC
plot(mlo.bAnom$VWC5.1, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$VWC5.1, mlo.bAnom$Overlap.postMed)

## Based on 10 day shallow SWC
plot(mlo.bAnom$VWC1.10, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$VWC1.10, mlo.bAnom$Overlap.postMed) ####
summary(lm(Overlap.postMed ~ VWC1.10, mlo.bAnom))

## Based on 10 day mid SWC
plot(mlo.bAnom$VWC3.10, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$VWC3.10, mlo.bAnom$Overlap.postMed)

## Based on 10 day deep SWC
plot(mlo.bAnom$VWC5.10, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$VWC5.10, mlo.bAnom$Overlap.postMed)

# Bouts with more diverse water use have higher LH flux ----
plot(mlo.bAnom$lh10, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$lh10, mlo.bAnom$Overlap.postMed)

# Bouts with more diverse water use have lower LH variation ----
plot(mlo.bAnom$lhVar, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$lhVar, mlo.bAnom$Overlap.postMed)

