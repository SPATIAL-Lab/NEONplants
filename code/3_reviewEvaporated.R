library(rgbif)

# Prep ####
## Read data
p = read.csv("data/plants.csv")
s = read.csv("data/soils.csv")

# Review Dex ####
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
abline(v = -30, lty = 3)

## Low values
p.low = p[p$Dex.off < -30,]

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
