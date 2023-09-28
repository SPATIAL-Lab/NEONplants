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
