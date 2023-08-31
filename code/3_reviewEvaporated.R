library(rgbif)

# Prep ####
## Read data
p = read.csv("data/plants.csv")

# Review Dex ####
## Calculate and distribution
p$Dex = p$d2H - p$d18O * 8
plot(density(p$Dex))

## Low values
p.low = p[p$Dex < -40,]

## Species and counts
s.low = data.frame(table(p.low$Species))

## Get common names from GBIF
s.low$Vernacular = rep("")
for(i in seq_along(s.low$Var1)){
  gbif = name_backbone(s.low$Var1[i])
  use = as.data.frame(name_usage(gbif$usageKey)$data)
  if(length(use$vernacularName > 0)){
    s.low$Vernacular[i] = use$vernacularName
  }
}

## Review - mostly shrubs and grasses
View(s.low)
