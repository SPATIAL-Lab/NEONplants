library(isoWater)

# Get soil and plant data from wiDB
d = wiDB_data(projects = '00384')
d = d$data

# Add bout index ----
## Extract bout index from sample ID
bout = substr(d$Sample_ID, 8, 12)
sort(unique(bout))

d$Bout = bout

## Remove one bout with no plant samples
d = d[d$Bout != "WOOD3",]

# Plant species ----
## Prep
p = d[d$Type == "Stem",]
p$Species = p$Sample_Comments

## Current comments
sort(unique(p$Species))

## Substitutions
p$Species = sub("Artemesia frigida", "Artemisia frigida", p$Species)
p$Species = sub("Artemisia tridentate", "Artemisia tridentata", p$Species)
p$Species = sub("Betula alleghensis", "Betula alleghaniensis", p$Species)
p$Species = sub("Bouteloua curtipendula \\(Michx.\\) Torr.", "Bouteloua curtipendula", p$Species)
p$Species = sub("Diosypros virginiana, all small trees/saplings", "Diosypros virginiana", p$Species)
p$Species = sub("Liquidambar styraciflua \\(alt.\\)", "Liquidambar styraciflua", p$Species)
p$Species = sub("Pinus palustris Mill.", "Pinus palustris", p$Species)
p$Species = sub("Pinus sabiniana Douglas ex Douglas", "Pinus sabiniana", p$Species)
p$Species = sub("Poaceae sp.; Achnatherum richardsonii", "Poaceae sp.", p$Species)
p$Species = sub("Poaceae sp.; Phleum pratense", "Poaceae sp.", p$Species)
p$Species = sub("Quercus douglasii Hook. & Arn.", "Quercus douglasii", p$Species)
p$Species = sub("Quercus falcata \\(alt.\\)", "Quercus falcata", p$Species)
p$Species = sub("Quercus stellata Wangenh.", "Quercus stellata", p$Species)
p$Species = sub("Quercus wislizeni A. DC.", "Quercus wislizeni", p$Species)
p$Species = sub("Schizachryium scoparium", "Schizachyrium scoparium", p$Species)
p$Species = sub("Smilax bona-nox L.", "Smilax bona-nox", p$Species)
p$Species = sub("Sorghastrum nutans \\(L.\\) Nash", "Sorghastrum nutans", p$Species)

## Check
sort(unique(p$Species))

# Soil depths ----
## Have a look
s = d[d$Type == "Soil",]
sort(unique(s$Depth_meters))

## Categorical
s$Depth = rep("")
s$Depth[s$Depth_meters <= 0.1] = "0-10 cm"
s$Depth[s$Depth_meters > 0.1 & s$Depth_meters <= 0.2] = "10-20 cm"
s$Depth[s$Depth_meters > 0.2 & s$Depth_meters <= 0.4] = "20-40 cm"
s$Depth[s$Depth_meters > 0.4] = "40+ cm"

## Check
sort(unique(s$Depth))

# Remove a few outliers ----
s = s[!(s$Bout == "CLBJ4" & s$d18O < -10),]
s = s[!(s$Bout == "CPER1" & s$d2H > -40),]
s = s[!(s$Bout == "SCBI6" & s$d2H < -140),]
s = s[!(s$Bout == "WOOD7" & s$d2H < -140),]

p = p[!(p$Bout == "CPER2" & p$d18O > 20),]
p = p[!(p$Bout == "WOOD6" & p$d18O < -15),]

# Save ----
write.csv(p, file = "data/plants.csv", row.names = FALSE)
write.csv(s, file = "data/soils.csv", row.names = FALSE)
