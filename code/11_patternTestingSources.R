library(terra)

## Water source info
ssd = read.csv("out/ssd.csv")

## Properties
bouts = read.csv("out/Bouts.csv")
sites = values(vect("out/Sites.shp"))

## Work with bouts w/ 10 or more individuals
ssd = ssd[ssd$n >= 10,]

# Wetter sites have less diverse water use ----
## Add site means to sites
ssd$Site = substr(ssd$Bout, 1, 4)
sites$SSD50 = rep(0)
for(i in seq_along(sites$Site)){
  sites$SSD50[i] = 
    mean(ssd$SSD50[ssd$Site == sites$Site[i]])
}

## Add all bouts to mlo.sites
ssd.sites = merge(ssd, sites)
boxplot(SSD50 ~ Site, ssd)

## Based on MAP - no
plot(sites$MAP, sites$SSD50)

## Based on VPD - yes
plot(sites$VPDmax, sites$SSD50)

# Wetter bouts have less diverse water use ----
## Add bouts to mlo.bouts
ssd.bouts = merge(ssd, bouts, by = "Bout")

## Based on 10 day precip
plot(ssd.bouts$Precip10, ssd.bouts$SSD50)

## Based on 30 day precip
plot(ssd.bouts$Precip30, ssd.bouts$SSD50)

## Based on 60 day precip
plot(ssd.bouts$Precip60, ssd.bouts$SSD50)

## Based on 1 day shallow SWC
plot(ssd.bouts$VWC1.1, ssd.bouts$SSD50)

## Based on 1 day mid SWC
plot(ssd.bouts$VWC3.1, ssd.bouts$SSD50)

## Based on 1 day deep SWC
plot(ssd.bouts$VWC5.1, ssd.bouts$SSD50)

## Based on 10 day shallow SWC
plot(ssd.bouts$VWC1.10, ssd.bouts$SSD50)

## Based on 10 day mid SWC
plot(ssd.bouts$VWC3.10, ssd.bouts$SSD50)

## Based on 10 day deep SWC
plot(ssd.bouts$VWC5.10, ssd.bouts$SSD50)

# Bouts with more diverse water use have higher LH flux ----
plot(ssd.bouts$lh10, ssd.bouts$SSD50)

# Bouts with more diverse water use have lower LH variation ----
plot(ssd.bouts$lhVar, ssd.bouts$SSD50)

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

