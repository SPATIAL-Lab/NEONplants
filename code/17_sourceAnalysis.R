library(terra)

# Summarize source depths for all individuals
## Load data
load("out/sourceStats.rda")

## Get all species and bouts
species = bouts = character()
for(i in 1:length(source.stats)){
  species = c(species, source.stats[[i]]$taxa)
  bouts = c(bouts, rep(names(source.stats)[i], length(source.stats[[i]]$taxa)))
}

## Extract median of mean depth
rd.50 = numeric()
for(i in 1:length(source.stats)){
  rdi.50 = apply(source.stats[[i]]$source.depths, 1, median)
  rd.50 = c(rd.50, rdi.50)
}

## Compile results
rd = data.frame("Bout" = bouts, "Species" = species, "Med" = rd.50)
plot(density(rd$Med))

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
rd.50 = rd.IQR = rd.n = numeric()
for(i in 1:length(sp)){
  rd.sub = rd$Med[rd$Species == sp[i]]
  rd.50 = c(rd.50, median(rd.sub))
  rd.IQR = c(rd.IQR, IQR(rd.sub))
  rd.n = c(rd.n, length(rd.sub))
}

## Compile results
rd.sp = data.frame("Species" = sp, "Form" = gf, "n" = rd.n, 
                   "Med" = rd.50, "IQR" = rd.IQR)
rd.sp.med = rd.sp[order(rd.sp$Med),]
rd.sp.iqr = rd.sp[order(rd.sp$IQR),]
rd.sp.iqr = rd.sp.iqr[rd.sp.iqr$n >= 5,]

## Have a look
head(rd.sp.med)
tail(rd.sp.med)
head(rd.sp.iqr)
tail(rd.sp.iqr)

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

barplot(iqr.plot$IQR, names = "", xlim = c(0.2, 12.2),
        main = "Source depth IQR", axes = FALSE, width = 1, 
        space = 0.1, col = "seagreen")
axis(2, las = 1)
for(i in c(1:5, 7:11)){
  text(0.6 + 1.1 * (i - 1), -0.05, bquote(italic(.(iqr.plot$Species[i]))), 
       srt = 45, xpd = NA, adj = 1)
}

dev.off()

# Analyze source depth by growth form
## Summarize by growth form
rd.50 = rd.25 = rd.75 = rd.n = numeric()
gf = unique(gf)
for(i in 1:length(gf)){
  rd.sub = rd.sp$Med[rd.sp$Form == gf[i]]
  rd.50 = c(rd.50, median(rd.sub))
  rd.25 = c(rd.25, quantile(rd.sub, 0.25))
  rd.75 = c(rd.75, quantile(rd.sub, 0.75))
  rd.n = c(rd.n, length(rd.sub))
}

## Compile results and look
rd.gf = data.frame("Form" = gf, "n" = rd.n, "Med" = rd.50, 
                   ".25" = rd.25, ".75" = rd.75)
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

## Have a look - nothing 
plot(bouts[, c(4:6, 15, 16)])
plot(bouts[, c(7:12, 15, 16)])
plot(bouts[, 13:16])

# Analyze source depth by bout properties - species level data
## Average median depth for species per bout
bouts$n.sp = bouts$mean.sp = bouts$sd.sp = numeric(nrow(bouts))
for(i in seq_along(bouts$Site)){
  rd.sub = rd[rd$Site == bouts$Site[i],]
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

## Have a look - mean depth increases and sd decreases as soil is wetter 
plot(bouts[, c(4:6, 17, 18)])
plot(bouts[, c(7:12, 17, 18)])
summary(lm(mean.sp ~ VWC3.10, bouts))
summary(lm(mean.sp ~ VWC1.10, bouts))
summary(lm(sd.sp ~ VWC3.10 , bouts))
summary(lm(sd.sp ~ VWC1.10 , bouts))
plot(bouts[, c(13:14, 17, 18)])
