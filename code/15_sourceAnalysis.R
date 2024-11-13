library(terra)

# Summarize source depths for all individuals
## Load data
load("out/sourceStats.rda")

## Get all species and bouts and KLds
species = bouts = character()
kld = numeric()
for(i in 1:length(source.stats)){
  species = c(species, source.stats[[i]]$taxa)
  bouts = c(bouts, rep(names(source.stats)[i], length(source.stats[[i]]$taxa)))
  kld = c(kld, source.stats[[i]]$KLd)
}

## Extract median of mean depth
rd.IQR = rd.50 = numeric()
for(i in 1:length(source.stats)){
  rdi.50 = apply(source.stats[[i]]$source.depths, 1, median)
  rd.50 = c(rd.50, rdi.50)
  rdi.IQR = apply(source.stats[[i]]$source.depths, 1, IQR)
  rd.IQR = c(rd.IQR, rdi.IQR)
}

## Compile results
rd = data.frame("Bout" = bouts, "Species" = species, "Med" = rd.50,
                "IQR" = rd.IQR, "KLd" = kld)

# Analyze source depth by species
## Unique species
sp = unique(species)

## Growth form...need trait data but limited coverage in TRY
gf = c("grass", "grass", 
       "shrub", "shrub", 
       "grass", "grass",
       "tree", "grass",
       "tree", "tree",
       "tree", "tree",
       "tree", "tree",
       "tree", "tree",
       "tree", "grass",
       "grass", "grass",
       "tree", "tree",
       "tree", "tree",
       "tree", "tree",
       "shrub", "grass",
       "grass", "shrub",
       "shrub", "tree",
       "tree", "grass",
       "tree", "tree",
       "tree", "tree",
       "tree", "shrub",
       "grass", "tree",
       "shrub", "tree",
       "tree", "tree",
       "grass", "grass",
       "grass", "shrub",
       "tree")

## Summarize by species
rd.50 = rd.IQR = rd.n = rd.IQR.mean = rd.KLd = numeric()
for(i in 1:length(sp)){
  rd.sub = rd[rd$Species == sp[i],]
  rd.50 = c(rd.50, median(rd.sub$Med))
  ### IQR of the median depths per species
  rd.IQR = c(rd.IQR, IQR(rd.sub$Med))
  ### Species mean IQR of individual source depths
  rd.IQR.mean = c(rd.IQR.mean, mean(rd.sub$IQR))
  rd.n = c(rd.n, length(rd.sub$Med))
  ### Species mean KLd
  rd.KLd = c(rd.KLd, mean(rd.sub$KLd))
}

## Compile results
rd.sp = data.frame("Species" = sp, "Form" = gf, "n" = rd.n, "KLd" = rd.KLd,
                   "Med" = rd.50, "IQR" = rd.IQR, "IQR.mean" = rd.IQR.mean)
rd.sp.med = rd.sp[order(rd.sp$Med),]
rd.sp.iqr = rd.sp[order(rd.sp$IQR),]
rd.sp.iqr = rd.sp.iqr[rd.sp.iqr$n >= 5,]

## Have a look
head(rd.sp.med)
tail(rd.sp.med)
head(rd.sp.iqr)
tail(rd.sp.iqr)

## Intermediate source depths are slightly less well constrained
plot(rd.sp$Med, rd.sp$IQR.mean)
plot(rd.sp$Med, rd.sp$KLd)

## Plot species highest/lowest
blank = rd.sp[1,]
blank[1,] = rep(NA)
med.plot = rbind(head(rd.sp.med, 5), blank, tail(rd.sp.med, 5))
iqr.plot = rbind(head(rd.sp.iqr, 5), blank, tail(rd.sp.iqr, 5))

png("out/sourceSpecies.png", width = 10, height = 5,
    units = "in", res = 600)
layout(matrix(c(1, 2), ncol = 2))

par(mar = c(8, 5, 3, 1))
barplot(med.plot$Med - 1, names = "", xlim = c(0.2, 12.2), 
        ylim = c(4, 0), axes = FALSE, main = "Median source depth",
        width = 1, space = 0.1, col = "seagreen")
axis(2, 4:0, c("GW", ">40 cm", "20-40 cm", "10-20 cm", "0-10 cm"), las = 1)
for(i in c(1:5, 7:11)){
  text(0.6 + 1.1 * (i - 1), 4.2, bquote(italic(.(med.plot$Species[i]))), 
       srt = 45, xpd = NA, adj = 1)
}
for(i in c(1:5, 7:11)){
  text(0.6 + 1.1 * (i - 1), 3.85, round(med.plot$KLd[i], 1), 
       srt = 90, adj = c(0, 0.5))
}


barplot(iqr.plot$IQR, names = "", xlim = c(0.2, 12.2),
        main = "Source depth IQR", axes = FALSE, width = 1, 
        space = 0.1, col = "seagreen")
axis(2, las = 1)
for(i in c(1:5, 7:11)){
  text(0.6 + 1.1 * (i - 1), -0.05, bquote(italic(.(iqr.plot$Species[i]))), 
       srt = 45, xpd = NA, adj = 1)
}
for(i in c(1:5, 7:11)){
  if(i == 1) {y = 0.1} else {y = 0.04}
  text(0.6 + 1.1 * (i - 1), y, round(iqr.plot$KLd[i], 1), 
       srt = 90, adj = c(0, 0.5))
}

dev.off()

# Analyze source depth by growth form
## Summarize by growth form
rd.50 = rd.25 = rd.75 = rd.IQR.mean = rd.n = rd.KLd = numeric()
gf = unique(gf)
for(i in 1:length(gf)){
  rd.sub = rd.sp[rd.sp$Form == gf[i],]
  rd.50 = c(rd.50, median(rd.sub$Med))
  rd.25 = c(rd.25, quantile(rd.sub$Med, 0.25))
  rd.75 = c(rd.75, quantile(rd.sub$Med, 0.75))
  rd.IQR.mean = c(rd.IQR.mean, mean(rd.sub$IQR.mean))
  rd.n = c(rd.n, length(rd.sub$Med))
  rd.KLd = c(rd.KLd, mean(rd.sub$KLd))
}

## Compile results and look
rd.gf = data.frame("Form" = gf, "n" = rd.n, "KLd" = rd.KLd, "Med" = rd.50, 
                   ".25" = rd.25, ".75" = rd.75, "IQR.mean" = rd.IQR.mean)
rd.gf

## Plot by GF
png("out/sourceGF.png", width = 6, height = 4, units = "in", res = 600)
par(mar = c(3, 8, 1, 1))
plot(rd.gf$Med, pch = 21, bg = "seagreen", cex = 3, axes = FALSE, 
     ylim = c(3.3, 2), xlim = c(0.8, 3.2), xlab = "", 
     ylab = "", lwd = 2)
arrows(1:3, rd.gf$X.25, y1 = rd.gf$X.75, length = 0, lwd = 2)
points(rd.gf$Med, pch = 21, bg = "seagreen", cex = 3, lwd = 2)
axis(1, 1:3, c("Grass", "Shrub", "Tree"))
axis(2, c(2, 3.3), labels = FALSE, lwd.ticks = 0)
axis(2, 2:3, c("10-20 cm", "20-40 cm"), las = 1)
mtext("Median source depth", 2, 5.5)
text(1:3, 2.1, round(rd.gf$KLd, 1), srt = 90, adj = c(1, 0.5))

dev.off()

# Analyze source depth by site properties - individual level data
## Properties
sites = values(vect("out/Sites.shp"))
rd$Site = substr(rd$Bout, 1, 4)

## Summarize individual depths by site
sites$mean = sites$sd = numeric(nrow(sites))
for(i in seq_along(sites$Site)){
  sites$mean[i] = mean(rd$Med[rd$Site == sites$Site[i]]) 
  sites$sd[i] = sd(rd$Med[rd$Site == sites$Site[i]]) 
}

## Have a look - mean depth marginally increases w/ MAP, sd increases w/ VPDmax
plot(sites)
summary(lm(mean ~ MAP, sites))
summary(lm(sd ~ VPDmax, sites))

# Analyze source depth by site properties - species level data
## Average median depth for species per site
sites$n.sp = sites$mean.sp = sites$sd.sp = numeric(nrow(sites))
for(i in seq_along(sites$Site)){
  rd.sub = rd[rd$Site == sites$Site[i],]
  sp.sub = unique(rd.sub$Species)
  n.sp = length(sp.sub)
  sp.mean = numeric(n.sp)
  for(j in seq_along(sp.sub)){
    sp.mean[j] = mean(rd.sub$Med[rd.sub$Species == sp.sub[j]])
  }
  sites$n.sp[i] = n.sp
  sites$mean.sp[i] = mean(sp.mean)
  if(n.sp > 1){
    sites$sd.sp[i] = sd(sp.mean)
  }
}

## Have a look - same responses, weaker
plot(sites)
summary(lm(mean.sp ~ MAP, sites))
summary(lm(sd.sp ~ VPDmax, sites))

# Analyze source depth by bout properties - individual level data
bouts = read.csv("out/Bouts.csv")

## Summarize rd per bout
bouts$mean = bouts$sd = numeric(nrow(bouts))
for(i in seq_along(bouts$Bout)){
  rd.sub = rd[rd$Bout == bouts$Bout[i],]
  bouts$mean[i] = mean(rd.sub$Med)
  bouts$sd[i] = sd(rd.sub$Med)
}

## Have a look - nothing much, maybe lower sd w/ wetter soil
plot(bouts[, c(4:6, 15, 16)])
plot(bouts[, c(7:12, 15, 16)])
summary(lm(sd ~ VWC3.10, bouts))
plot(bouts[, 13:16])

# Analyze source depth by bout properties - species level data
## Average median depth for species per bout
bouts$n.sp = bouts$mean.sp = bouts$sd.sp = numeric(nrow(bouts))
for(i in seq_along(bouts$Site)){
  rd.sub = rd[rd$Bout == bouts$Bout[i],]
  sp.sub = unique(rd.sub$Species)
  n.sp = length(sp.sub)
  sp.mean = numeric(n.sp)
  for(j in seq_along(sp.sub)){
    sp.mean[j] = mean(rd.sub$Med[rd.sub$Species == sp.sub[j]])
  }
  bouts$n.sp[i] = n.sp
  bouts$mean.sp[i] = mean(sp.mean)
  if(n.sp > 1){
    bouts$sd.sp[i] = sd(sp.mean)
  }
}

## Have a look 
plot(bouts[, c(4:6, 17, 18)])
plot(bouts[, c(7:12, 17, 18)])
plot(bouts[, c(13:14, 17, 18)])

# Are individual or species mean depth differences w/in bouts related to 
# source range? Not really.
ssd.stats = read.csv("out/ssd.csv")
bouts$srangeO = ssd.stats$srangeO[match(bouts$Bout, ssd.stats$Bout)]
bouts$SSD50 = ssd.stats$SSD50[match(bouts$Bout, ssd.stats$Bout)]

summary(lm(sd ~ srangeO, bouts[bouts$srangeO > 4,]))
summary(lm(sd.sp ~ srangeO, bouts[bouts$srangeO > 4,]))

bouts$sd.norm = bouts$sd / bouts$SSD50
plot(sd.norm ~ srangeO, bouts)
summary(lm(sd.norm ~ srangeO, bouts))
bouts$sd.norm = bouts$sd / bouts$SSD50

png("out/sdXsw.png", width = 9, height = 5, units = "in", res = 600)
layout(matrix(c(1, 2), ncol = 2))
plot(bouts$srangeO, bouts$sd, pch = 21, cex = 2, bg = "seagreen",
     xlab = expression("Source water "*delta*{18}*"O range"),
     ylab = "SD of median source depth (individuals)")
plot(bouts$srangeO, bouts$sd.sp, pch = 21, cex = 2, bg = "seagreen",
     xlab = expression("Source water "*delta*{18}*"O range"),
     ylab = "SD of median source depth (species)")

dev.off()

# Plot site level responses
png("out/sourceSites.png", width = 9, height = 5,
    units = "in", res = 600)
layout(matrix(c(1, 2), ncol = 2))

par(mar = c(5, 5, 3, 1))
plot(sites$MAP, sites$mean, ylim = c(3.7, 2), xlab = "Site MAP (mm)",
     ylab = "", main = "Mean source depth", axes = FALSE)
abline(lm(mean ~ MAP, sites), lty = 2)
points(sites$MAP, sites$mean, pch = 21, cex = 2, bg = "seagreen",
       lwd = 2)
axis(2, c(2, 3.7), labels = FALSE, lwd.ticks = 0)
axis(2, 2:3, c("10-20 cm", "20-40 cm"), las = 1)
axis(1)
x = par("usr")[2] - 0.05 * diff(range(par("usr")[1:2]))
y = par("usr")[4] + 0.1 * diff(range(par("usr")[3:4]))
text(x, y, "p = 0.06", pos = 2)

plot(sites$VPDmax, sites$sd, xlab = "Site VPDmax (hPa)", ylab = "", 
     ylim = c(0.25, 0.95), main = "Source depth SD", axes = FALSE)
abline(lm(sd ~ VPDmax, sites), lty = 2)
points(sites$VPDmax, sites$sd, pch = 21, cex = 2, bg = "seagreen",
       lwd = 2)
axis(2, c(0.3, 0.6, 0.9), las = 1)
axis(1)
x = par("usr")[2] - 0.05 * diff(range(par("usr")[1:2]))
y = par("usr")[3] + 0.1 * diff(range(par("usr")[3:4]))
text(x, y, "p = 0.04", pos = 2)

dev.off()
