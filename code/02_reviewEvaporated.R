library(rgbif)

# Prep ----
## Read data
p = read.csv("data/plants.csv")
s = read.csv("data/soils.csv")
irms = read.csv("data/irms.csv")

# Review Dex ----
## Calculate and distribution
p$Dex = p$d2H - p$d18O * 8
plot(density(p$Dex))

## Minimum Dex soil values per bout
bouts = unique(s$Bout)
s.dex.min = numeric(length(bouts))
for(i in seq_along(bouts)){
  s.bout = s[s$Bout == bouts[i],]
  s.dex.min[i] = min(s.bout$d2H - 8 * s.bout$d18O)
}

## Compare plant values with minimum soil values
p$Dex.off = p$Dex - s.dex.min[match(p$Bout, bouts)]
plot(density(p$Dex.off))
low = -30
abline(v = low, lty = 3)

## Low values
p.low = p[p$Dex.off < low,]

## Species and counts
sp.low = data.frame(table(p.low$Species))

## Get common names from GBIF
sp.low$Vernacular = rep("")
for(i in seq_along(sp.low$Var1)){
  gbif = name_backbone(sp.low$Var1[i])
  use = as.data.frame(name_usage(gbif$usageKey)$data)
  if(length(use$vernacularName > 0)){
    sp.low$Vernacular[i] = use$vernacularName
  }
}

## Review - mostly shrubs and grasses
View(sp.low)

## Save
write.csv(p, "out/plantWDex.csv", row.names = FALSE)

# Compare with irms data ----
s$ID = substr(s$Sample_ID, 8, nchar(s$Sample_ID))
p$ID = substr(p$Sample_ID, 8, nchar(p$Sample_ID))
sirms = merge(s, irms, by = "ID")
pirms = merge(p, irms, by = "ID")

## Soil plots
plot(sirms$d2H.irms, sirms$d2H)
abline(0, 1)
plot(sirms$d18O.irms, sirms$d18O)
abline(0, 1)

## Soil stats
s.diff.h = apply(cbind(sirms$d2H, sirms$d2H.irms), 1, diff)
s.diff.o = apply(cbind(sirms$d18O, sirms$d18O.irms), 1, diff)
s.diff.o = s.diff.o[s.diff.o < 5]
t.test(s.diff.h)
t.test(s.diff.o)
sd(s.diff.h)
sd(s.diff.o)

## Plant plots
plot(pirms$d2H.irms, pirms$d2H)
abline(0, 1)
plot(pirms$d18O.irms, pirms$d18O)
abline(0, 1)

## Plant stats
p.diff.h = apply(cbind(pirms$d2H, pirms$d2H.irms), 1, diff)
p.diff.o = apply(cbind(pirms$d18O, pirms$d18O.irms), 1, diff)
t.test(p.diff.h)
t.test(p.diff.o)
sd(p.diff.h)
sd(p.diff.o)
plot(density(p.diff.o))

## View the plant results for O
View(data.frame(pirms$ID, pirms$Species, p.diff.o))

## Write out these results
pirms = cbind(pirms, p.diff.h, p.diff.o)
write.csv(pirms, "out/plantDiff.csv")

## Plot all
wh = 8
widths = sapply(wh * 2.54 * c(0.08, 0.46, 0.46), lcm)
png("out/IRMS.png", wh, wh, units = "in", res = 600)
layout(matrix(1:9, nrow = 3), widths = widths,
       heights = widths)

par(mar = rep(0, 4))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")

par(mar = c(5, 0, 1, 0))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
text(0, 0, expression(bold("Soil")), cex = 1.75, srt = 90, 
     col = "salmon2")

plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
text(0, 0, expression(bold("Xylem")), cex = 1.75, srt = 90, 
     col = "seagreen")

par(mar = c(0, 5, 0, 1))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
text(0, 0, expression(bold(delta^2*"H (VSMOW)")), cex = 1.75)

par(mar = c(5, 5, 1, 1))
plot(sirms$d2H.irms, sirms$d2H, pch = 21, bg = "salmon2",
     xlab = "IRMS", ylab = "CRDS", cex.lab = 1.5)
abline(0, 1, lwd = 2)
points(sirms$d2H.irms, sirms$d2H, pch = 21, bg = "salmon2", cex = 1.5)
text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
     paste0("RMSE = ", round(sqrt(mean(s.diff.h ^ 2)), 1), "\u2030"), 
     adj = c(1, 0))

plot(pirms$d2H.irms, pirms$d2H, pch = 21, bg = "seagreen",
     xlab = "IRMS", ylab = "CRDS", cex.lab = 1.5)
abline(0, 1)
points(pirms$d2H.irms, pirms$d2H, pch = 21, bg = "seagreen", cex = 1.5)
text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
     paste0("RMSE = ", round(sqrt(mean(p.diff.h ^ 2)), 1), "\u2030"), 
     adj = c(1, 0))

par(mar = c(0, 5, 0, 1))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
text(0, 0, expression(bold(delta^{18}*"O (VSMOW)")), cex = 1.75)

par(mar = c(5, 5, 1, 1))
plot(sirms$d18O.irms, sirms$d18O, pch = 21, bg = "salmon2", xlim = c(-9.5, -2),
     xlab = "IRMS", ylab = "CRDS", cex.lab = 1.5)
abline(0, 1)
points(sirms$d18O.irms, sirms$d18O, pch = 21, bg = "salmon2", cex = 1.5)
text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
     paste0("RMSE = ", round(sqrt(mean(s.diff.o ^ 2)), 1), "\u2030"), 
     adj = c(1, 0))

plot(pirms$d18O.irms, pirms$d18O, pch = 21, bg = "seagreen",
     xlab = "IRMS", ylab = "CRDS", cex.lab = 1.5)
abline(0, 1)
points(pirms$d18O.irms, pirms$d18O, pch = 21, bg = "seagreen", cex = 1.5)
text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
     paste0("RMSE = ", round(sqrt(mean(p.diff.o ^ 2)), 1), "\u2030"), 
     adj = c(1, 0))

dev.off()

