library(rgbif)

# Pull raw data files -not portable- ----
## Read data
p = read.csv("out/plantWDex.csv")
p.diff = read.csv("out/plantDiff.csv")

## Get list of job numbers
sid = strsplit(p$Sample_ID, "_")
sid = sapply(sid, "[[", 1)
sid = unique(sid)

## Directories 
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

# Process metrics -not portable- ----
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

write.csv(sshift, "out/sshift.csv")

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

## Fraction of data affected
sum(is.na(p$SlopeShift))
sum(p$SlopeShift < -15) / nrow(p)

## Screened data w/ IRMS offsets
p.diff.ok = p.diff[p.diff$SlopeShift > -15,]

## Are averages for the screened data different than zero?
t.test(p.diff.ok$d2H - p.diff.ok$d2H.irms)
t.test(p.diff.ok$d18O - p.diff.ok$d18O.irms)

## Offset correction for d18O
p.diff.ok$d18O.oc = p.diff.ok$d18O - mean(p.diff.ok$d18O - p.diff.ok$d18O.irms)
p.diff$d18O.oc = p.diff$d18O - mean(p.diff.ok$d18O - p.diff.ok$d18O.irms)
p$d18O.oc = p$d18O - mean(p.diff.ok$d18O - p.diff.ok$d18O.irms)

png("out/Dex_screening.png")
plot(density(p$Dex.off))
lines(density(p$Dex.off[p$SlopeShift > -15]), col = "blue")
lines(density(p$Dex.off[p$SlopeShift > -5]), col = "red")
dev.off()

## Plot p.diff for normal and low slopeshifts
png("out/slopeShiftIRMS.png", 8, 4, units = "in", res = 600)
layout(matrix(1:2, nrow = 1))
par(mar = c(5, 5, 1, 1))

plot(p.diff$d2H.irms, p.diff$d2H, pch = 21, bg = "gray90",
     xlab = expression(delta^2*"H IRMS"), 
     ylab = expression(delta^2*"H CRDS"))
abline(0, 1)
points(p.diff$d2H.irms, p.diff$d2H, pch = 21, bg = "gray90")
points(p.diff.ok$d2H.irms, p.diff.ok$d2H, 
       pch = 21, bg = "seagreen")
text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
     paste0("RMSE = ", 
            round(sqrt(mean((p.diff.ok$d2H - p.diff.ok$d2H.irms) ^ 2)), 1), 
            "\u2030"), adj = c(1, 0))

plot(p.diff$d18O.irms, p.diff$d18O, pch = 21, bg = "gray90",
     xlab = expression(delta^{18}*"O IRMS"), 
     ylab = expression(delta^{18}*"O CRDS"),
     ylim = range(c(p.diff$d18O, p.diff$d18O.oc)))
abline(0, 1)
points(p.diff$d18O.irms, p.diff$d18O, pch = 21, bg = "gray90")
points(p.diff.ok$d18O.irms, p.diff.ok$d18O.oc, 
       pch = 21, bg = "seagreen")
text(par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
     par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
     paste0("RMSE = ", 
            round(sqrt(mean((p.diff.ok$d18O.oc - p.diff.ok$d18O.irms) ^ 2)), 1), 
            "\u2030"), adj = c(1, 0))

dev.off()

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
specSpec = specSpec[rev(order(specSpec$LowFrac)),]
write.csv(specSpec, "out/screenedBySpecies.csv")

## Plot highest and lowest
blank = specSpec[1,]
blank[1,] = rep(NA)
specplot = rbind(head(specSpec, 5), blank, tail(specSpec, 5))

png("out/specShiftSpecies.png", 6, 4, units = "in", res = 600)
par(mar = c(8, 5, 1, 1))
barplot(specplot$LowFrac, names = "", xlim = c(0.2, 12.2), 
        width = 1, space = 0.1, col = "seagreen")
mtext("Low SS Fraction", 2, 3)
for(i in c(1:5, 7:11)){
  text(0.6 + 1.1 * (i - 1), -0.1, bquote(italic(.(specplot$Species[i]))), 
       srt = 45, xpd = NA, adj = 1)
}

dev.off()

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
