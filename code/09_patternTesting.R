library(terra)

## Water source info
mlo = read.csv("out/MLO.csv")
mix = read.csv("out/mixStats.csv")

## Properties
bouts = read.csv("out/Bouts.csv")
sites = values(vect("out/Sites.shp"))

# Wetter sites have less diverse water use ----
## Set up
mlo$Site = substr(mlo$Bout, 1, 4)
mlo.sites = merge(mlo, sites)
boxplot(Overlap.raw ~ Site, mlo)
boxplot(Overlap.post ~ Site, mlo)
boxplot(Overlap.postMed ~ Site, mlo)

## Based on MAP - no
plot(mlo.sites$MAP, mlo.sites$Overlap.postMed)
summary(lm(Overlap.postMed ~ MAP, mlo.sites))

## Based on VPD - no
plot(mlo.sites$VPDmax, mlo.sites$Overlap.postMed)
summary(lm(Overlap.postMed ~ VPDmax, mlo.sites))

# Wetter bouts have less diverse water use ----
## Set up
mlo.bouts = merge(mlo, bouts, by = "Bout")

## Based on 10 day precip
plot(mlo.bouts$Precip10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$Precip10, mlo.bouts$Overlap.post)
plot(mlo.bouts$Precip10, mlo.bouts$Overlap.postMed)

## Based on 30 day precip
plot(mlo.bouts$Precip30, mlo.bouts$Overlap.raw)
plot(mlo.bouts$Precip30, mlo.bouts$Overlap.post)
plot(mlo.bouts$Precip30, mlo.bouts$Overlap.postMed)

## Based on 60 day precip
plot(mlo.bouts$Precip60, mlo.bouts$Overlap.raw)
plot(mlo.bouts$Precip60, mlo.bouts$Overlap.post)
plot(mlo.bouts$Precip60, mlo.bouts$Overlap.postMed)

## Based on 1 day shallow SWC
plot(mlo.bouts$VWC1.1, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC1.1, mlo.bouts$Overlap.post)
summary(lm(Overlap.post ~ VWC1.1, mlo.bouts))
plot(mlo.bouts$VWC1.1, mlo.bouts$Overlap.postMed)

## Based on 1 day mid SWC
plot(mlo.bouts$VWC3.1, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC3.1, mlo.bouts$Overlap.post)
summary(lm(Overlap.post ~ VWC3.1, mlo.bouts))
plot(mlo.bouts$VWC3.1, mlo.bouts$Overlap.postMed)

## Based on 1 day deep SWC
plot(mlo.bouts$VWC5.1, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC5.1, mlo.bouts$Overlap.post)
summary(lm(Overlap.post ~ VWC5.1, mlo.bouts))
plot(mlo.bouts$VWC5.1, mlo.bouts$Overlap.postMed)

## Based on 10 day shallow SWC
plot(mlo.bouts$VWC1.10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC1.10, mlo.bouts$Overlap.post)
summary(lm(Overlap.post ~ VWC1.10, mlo.bouts))
plot(mlo.bouts$VWC1.10, mlo.bouts$Overlap.postMed)

## Based on 10 day mid SWC
plot(mlo.bouts$VWC3.10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC3.10, mlo.bouts$Overlap.post)
summary(lm(Overlap.post ~ VWC3.10, mlo.bouts))
plot(mlo.bouts$VWC3.10, mlo.bouts$Overlap.postMed)
summary(lm(Overlap.postMed ~ VWC3.10, mlo.bouts))

## Based on 10 day deep SWC
plot(mlo.bouts$VWC5.10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$VWC5.10, mlo.bouts$Overlap.post)
summary(lm(Overlap.post ~ VWC5.10, mlo.bouts))
plot(mlo.bouts$VWC5.10, mlo.bouts$Overlap.postMed)

# Bouts with more diverse water use have higher LH flux ----
plot(mlo.bouts$lh10, mlo.bouts$Overlap.raw)
plot(mlo.bouts$lh10, mlo.bouts$Overlap.post)
plot(mlo.bouts$lh10, mlo.bouts$Overlap.postMed)

# Bouts with more diverse water use have lower LH variation ----
plot(mlo.bouts$lhVar, mlo.bouts$Overlap.raw)
plot(mlo.bouts$lhVar, mlo.bouts$Overlap.post)
plot(mlo.bouts$lhVar, mlo.bouts$Overlap.postMed)

