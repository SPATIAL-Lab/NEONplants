library(rgbif)

# Pull raw data files ----
## Read data
p = read.csv("out/plantWDex.csv")
p.diff = read.csv("out/plantDiff.csv")

## Get list of job numbers
sid = strsplit(p$Sample_ID, "_")
sid = sapply(sid, "[[", 1)
sid = unique(sid)

## Directories -not portable-
wd = "~/../Dropbox"

## Find runfiles files from jobnumbers
cf = rf = character()
for(i in 1:length(sid)){
  rf = append(rf, list.files(file.path(wd, "hids2046", "runfiles"),
                             sid[i], full.names = TRUE))
  rf = append(rf, list.files(file.path(wd, "hids2052", "runfiles"), 
                             sid[i], full.names = TRUE))
}

## Parse out the directory and runfile name
fstring = strsplit(rf, "/")
fpath = sapply(fstring, "[[", 7)
fname = sapply(fstring, "[[", 9)

## From runfile name get date string
fdate = substr(fname, 1, 6)

## Find all coordinator files
for(i in seq_along(fdate)){
  cfs = list.files(file.path(wd, fpath[i]), fdate[i], full.names = TRUE)
  cf = append(cf, cfs[length(cfs)])
}

# Process metrics ----
sshift = data.frame("Vial" = numeric(), "SlopeShift" = numeric(),
                    "ID" = character())
for(i in seq_along(rf)){
  ## Read files
  ids = read.csv(rf[i])
  d = read.csv(cf[i])
  
  ## Get vial #s from Port
  d.vials = strsplit(d$Port, "-")
  d.vials = as.numeric(sapply(d.vials, "[[", 2))
  
  ## Average the slope shift per vial
  if(length(unique(d.vials)) > 4){
    spec = data.frame("Vial" = unique(d.vials), "SlopeShift" = rep(0))
    for(j in seq_along(spec$Vial)){
      spec$SlopeShift[j] = mean(d$slope_shift[d.vials == spec$Vial[j]],
                                na.rm = TRUE)
    }
    
    ## Standardize to the vial 4 average
    spec$SlopeShift = spec$SlopeShift - spec$SlopeShift[spec$Vial == 4]
    
    ## Add sample IDs
    ids$ID = paste(ids$Identifier.2, ids$Description, sep = "_")
    spec$ID = ids$ID[match(spec$Vial, ids$Vial)]
    
    ## Bind
    sshift = rbind(sshift, spec)
  }
}

# Condense and match ----
## Average per ID
ssAve = data.frame("ID" = unique(sshift$ID), "SlopeShift" = rep(0))
for(i in seq_along(ssAve$ID)){
  ssAve$SlopeShift[i] = mean(sshift$SlopeShift[sshift$ID == ssAve$ID[i]],
                             na.rm = TRUE)
}

## Plot
ssAve = ssAve[!is.na(ssAve$SlopeShift),]
plot(density(ssAve$SlopeShift))

## Add back to datasets
p$SlopeShift = ssAve$SlopeShift[match(p$Sample_ID, ssAve$ID)]
p.diff$SlopeShift = ssAve$SlopeShift[match(p.diff$Sample_ID, ssAve$ID)]

# Explore ----
plot(p.diff$SlopeShift, p.diff$p.diff.o)
plot(p.diff$SlopeShift, p.diff$p.diff.h)

plot(p$SlopeShift, p$Dex.off)

sum(is.na(p$SlopeShift))
sum(p$SlopeShift < -15, na.rm = TRUE)

plot(density(p$Dex.off))
lines(density(p$Dex.off[p$SlopeShift > -15]), col = "blue")
lines(density(p$Dex.off[p$SlopeShift > -5]), col = "red")

# Species summaries ---
## Space
specSpec = data.frame("Species" = unique(p$Species), "Count" = rep(0),
                      "LowSS" = rep(0))

## How many low SS?
for(i in seq_along(specSpec$Species)){
  psub = p[p$Species == specSpec$Species[i] & !is.na(p$SlopeShift),]
  specSpec$Count[i] = nrow(psub)
  specSpec$LowSS[i] = sum(psub$SlopeShift < -15)
}
specSpec$LowFrac = specSpec$LowSS / specSpec$Count

## Get common names from GBIF
specSpec$Vernacular = rep("")
for(i in seq_along(specSpec$Species)){
  gbif = name_backbone(specSpec$Species[i])
  use = as.data.frame(name_usage(gbif$usageKey)$data)
  if(length(use$vernacularName > 0)){
    specSpec$Vernacular[i] = use$vernacularName
  }
}

# Site summaries ----
## Space
siteSpec = data.frame("Site" = unique(p$Site_ID), "Count" = rep(0),
                      "LowSS" = rep(0))

## How many low SS?
for(i in seq_along(siteSpec$Site)){
  psub = p[p$Site_ID == siteSpec$Site[i] & !is.na(p$SlopeShift),]
  siteSpec$Count[i] = nrow(psub)
  siteSpec$LowSS[i] = sum(psub$SlopeShift < -15)
}
siteSpec$LowFrac = siteSpec$LowSS / siteSpec$Count

# Save ----
p = p[p$SlopeShift > -15,]
write.csv(p[p$SlopeShift > -15,], "out/plantsNoLow.csv", row.names = FALSE)
