library(terra)

## Water source info
mlo = read.csv("out/MLO.csv")

## Properties
bouts = read.csv("out/Bouts.csv")
sites = values(vect("out/Sites.shp"))

## Work with bouts w/ 3 or more taxa
mlo = mlo[mlo$SpeciesNum > 2,]
mlo$Site = substr(mlo$Bout, 1, 4)

## Plot by site
boxplot(Overlap.raw ~ Site, mlo)
boxplot(Overlap.postMed ~ Site, mlo)

# Wetter sites have less diverse water use ----
## Add site means to sites
sites$Overlap.raw = sites$Overlap.post = sites$Overlap.postMed = rep(0)
for(i in seq_along(sites$Site)){
  sites$Overlap.raw[i] = 
    mean(mlo$Overlap.raw[mlo$Site == sites$Site[i]])
  sites$Overlap.post[i] = 
    mean(mlo$Overlap.post[mlo$Site == sites$Site[i]])
  sites$Overlap.postMed[i] = 
    mean(mlo$Overlap.postMed[mlo$Site == sites$Site[i]])
}

## Remove SJER which lacks MLO values
sites = sites[sites$Site != "SJER",]

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
plot(mlo.bouts$VWC5.1, mlo.bouts$Overlap.postMed)

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
summary(lm(Overlap.postMed ~ VWC3.10, mlo.bAnom))

## Based on 10 day deep SWC
plot(mlo.bAnom$VWC5.10, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$VWC5.10, mlo.bAnom$Overlap.postMed)

# Bouts with more diverse water use have higher LH flux ----
plot(mlo.bAnom$lh10, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$lh10, mlo.bAnom$Overlap.postMed)

# Bouts with more diverse water use have lower LH variation ----
plot(mlo.bAnom$lhVar, mlo.bAnom$Overlap.raw)
plot(mlo.bAnom$lhVar, mlo.bAnom$Overlap.postMed)

# Plots
## Sites w/ higher T demand have more diverse water use
png("out/overlapVPD.png", width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
plot(sites$VPDmax, sites$Overlap.postMed, xlab = "Site VPDmax (hPa)",
     ylab = "Source overlap") 
mod = lm(Overlap.postMed ~ VPDmax, sites)
abline(mod, lty = 2)
points(sites$VPDmax, sites$Overlap.postMed, pch = 21, cex = 2,
       bg = "seagreen", lwd = 2)
p = signif(summary(mod)$coefficients[2, 4], 1)
text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
     par("usr")[4] - 0.1 * diff(par("usr")[3:4]),
     paste("p =", p), pos = 2)

dev.off()

## Anomalously wet bouts have less diverse water use - 60 day precip
png("out/DoverlapPrecip.png", width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
plot(mlo.bAnom$Precip60, mlo.bAnom$Overlap.postMed, 
     xlab = "Bout 60-day precipitation anomaly (mm/day)", 
     ylab = "Bout source overlap anomaly") 
mod = lm(Overlap.postMed ~ Precip60, mlo.bAnom)
abline(mod, lty = 2)
points(mlo.bAnom$Precip60, mlo.bAnom$Overlap.postMed, pch = 21, cex = 2,
       bg = "seagreen", lwd = 2)
p = signif(summary(mod)$coefficients[2, 4], 1)
text(par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.1 * diff(par("usr")[3:4]),
     paste("p =", p), pos = 4)

dev.off()

## Anomalously wet bouts have less diverse water use - 10 day shallow SWC
png("out/DoverlapSWC.png", width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
plot(mlo.bAnom$VWC1.10, mlo.bAnom$Overlap.postMed, 
     xlab = "Bout 10-day shallow SWC anomaly (%)", 
     ylab = "Bout source overlap anomaly") 
mod = lm(Overlap.postMed ~ VWC1.10, mlo.bAnom)
abline(mod, lty = 2)
points(mlo.bAnom$VWC1.10, mlo.bAnom$Overlap.postMed, pch = 21, cex = 2,
       bg = "seagreen", lwd = 2)
p = signif(summary(mod)$coefficients[2, 4], 1)
text(par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.1 * diff(par("usr")[3:4]),
     paste("p =", p), pos = 4)

dev.off()

png("out/overlapSWC.png", width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
plot(mlo.bouts$VWC3.1, mlo.bouts$Overlap.postMed, pch = 21, cex = 2,
     bg = "seagreen", lwd = 2,
     xlab = "Bout 1-day intermediate SWC (%)", 
     ylab = "Bout source overlap") 

dev.off()
