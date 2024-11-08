library(neonUtilities)
library(rhdf5)

# Prep ----
p = read.csv("out/plantsNoLow.csv")
bouts = unique(p$Bout)
dts = p$Collection_Date[match(bouts, p$Bout)]
sts = substr(bouts, 1, 4)
bouts = data.frame("Bout" = bouts, "Date" = as.POSIXct(dts), "Site" = sts)

# Precipitation ----
bouts$Precip60 = bouts$Precip30 = bouts$Precip10 = rep(0)

## Cycle through bouts
for(i in 1:nrow(bouts)){
  pre = loadByProduct("DP1.00006.001", bouts$Site[i], 
                      as.character(bouts$Date[i] - 60 * 3600 * 24), 
                      as.character(bouts$Date[i]), timeIndex = 30, 
                      check.size = FALSE)

  ## Precipitation across 10, 30, and 60 days prior to sampling
  pre = pre$PRIPRE_30min
  pre = pre[pre$startDateTime < bouts$Date[i],]
  pre60 = pre[pre$startDateTime > bouts$Date[i] - 60 * 3600 * 24,]
  pre30 = pre[pre$startDateTime > bouts$Date[i] - 30 * 3600 * 24,]
  pre10 = pre[pre$startDateTime > bouts$Date[i] - 10 * 3600 * 24,]

  ## Plot and sum 60 day
  if(sum(!(is.na(pre60$priPrecipBulk))) > 0){
    plot(as.POSIXct(pre60$startDateTime), pre60$priPrecipBulk,
         type = "l", main = bouts$Bout[i])
    bouts$Precip60[i] = sum(pre60$priPrecipBulk, na.rm = TRUE)
    ## Rescale to rate (mm per day), accounting for missing obs
    bouts$Precip60[i] = bouts$Precip60[i] / sum(!(is.na(pre60$priPrecipBulk))) * 48
  }else{
    bouts$Precip60[i] = NA
  }

  ## Sum 30 day
  if(sum(!(is.na(pre30$priPrecipBulk))) > 0){
    bouts$Precip30[i] = sum(pre30$priPrecipBulk, na.rm = TRUE)
    bouts$Precip30[i] = bouts$Precip30[i] / sum(!(is.na(pre30$priPrecipBulk))) * 48
  }else{
    bouts$Precip30[i] = NA
  }
  
  ## Sum 10 day
  if(sum(!(is.na(pre10$priPrecipBulk))) > 0){
    bouts$Precip10[i] = sum(pre10$priPrecipBulk, na.rm = TRUE)
    bouts$Precip10[i] = bouts$Precip10[i] / sum(!(is.na(pre10$priPrecipBulk))) * 48
  }else{
    bouts$Precip10[i] = NA
  }
}

# Soil moisture ----
bouts$VWC1.10 = bouts$VWC3.10 = bouts$VWC5.10 = rep(0)
bouts$VWC1.1 = bouts$VWC3.1 = bouts$VWC5.1 = rep(0)

## Cycle through bouts
for(i in 1:nrow(bouts)){
  vwc = loadByProduct("DP1.00094.001", bouts$Site[i], 
                      as.character(bouts$Date[i] - 30 * 3600 * 24), 
                      as.character(bouts$Date[i]), timeIndex = 30,
                      check.size = FALSE)
  vwc = vwc$SWS_30_minute
  
  ## Shallow, mid and deep SWC for preceding 10 days
  vwc = vwc[vwc$startDateTime > bouts$Date[i] - 10 * 3600 * 24,]
  vwc = vwc[vwc$startDateTime < bouts$Date[i],]
  vwc = vwc[vwc$VSWCFinalQF == 0,]
  vwc1 = vwc[vwc$verticalPosition == "501",]
  vwc3 = vwc[vwc$verticalPosition == "503",]
  vwc5 = vwc[vwc$verticalPosition == "505",]
  
  if(sum(!(is.na(vwc1$VSWCMean))) > 0 & sum(!(is.na(vwc3$VSWCMean))) > 0
     & sum(!(is.na(vwc5$VSWCMean))) > 0){
    plot(vwc1$startDateTime, vwc1$VSWCMean, main = bouts$Bout[i])
    points(vwc3$startDateTime, vwc3$VSWCMean, col = "red")
    points(vwc5$startDateTime, vwc5$VSWCMean, col = "blue")
    
    bouts$VWC1.10[i] = mean(vwc1$VSWCMean, na.rm = TRUE)
    bouts$VWC3.10[i] = mean(vwc3$VSWCMean, na.rm = TRUE)
    bouts$VWC5.10[i] = mean(vwc5$VSWCMean, na.rm = TRUE)
  }else{
    bouts$VWC1.10[i] = bouts$VWC3.10[i] = bouts$VWC5.10[i] = NA
  }
  
  ## Now for preceding 1 day
  vwc = vwc[vwc$startDateTime > bouts$Date[i] - 3600 * 24,]
  vwc1 = vwc[vwc$verticalPosition == "501",]
  vwc3 = vwc[vwc$verticalPosition == "503",]
  vwc5 = vwc[vwc$verticalPosition == "505",]
  
  if(sum(!(is.na(vwc1$VSWCMean))) > 0 & sum(!(is.na(vwc3$VSWCMean))) > 0
     & sum(!(is.na(vwc5$VSWCMean))) > 0){
    bouts$VWC1.1[i] = mean(vwc1$VSWCMean, na.rm = TRUE)
    bouts$VWC3.1[i] = mean(vwc3$VSWCMean, na.rm = TRUE)
    bouts$VWC5.1[i] = mean(vwc5$VSWCMean, na.rm = TRUE)
  }else{
    bouts$VWC1.1[i] = bouts$VWC3.1[i] = bouts$VWC5.1[i] = NA
  }
}

# Latent heat flux ----
bouts$lh10 = bouts$lhVar = rep(0)

## Cycle through bouts
for(i in 1:nrow(bouts)){
  ## Space for downloads
  dld = file.path(tempdir(), bouts$Bout[i])
  if(!dir.exists(dld)){
    dir.create(dld)
    ## Get EC bundle
    zipsByProduct("DP4.00200.001", bouts$Site[i], 
                  as.character(bouts$Date[i] - 60 * 3600 * 24),
                  as.character(bouts$Date[i]), check.size = FALSE,
                  timeIndex = 30, savepath = dld)
  }
  
  lh = stackEddy(file.path(dld, "filesToStack00200"))
  lh = lh[[1]]

  ## LH flux for preceding 60 days
  lh = lh[lh$timeBgn > bouts$Date[i] - 60 * 3600 * 24,]
  lh = lh[lh$timeBgn < bouts$Date[i],]
  lh = lh[lh$qfqm.fluxH2o.nsae.qfFinl == 0,]

  if(nrow(lh) > 0){
    plot(lh$timeBgn, lh$data.fluxH2o.nsae.flux, main = bouts$Bout[i])
    
    ## Get max values for each day
    days = unique(as.Date(lh$timeBgn))
    days = data.frame("Day" = days, "LH.max" = rep(0))
    for(j in seq_along(days$Day)){
      day = lh$data.fluxH2o.nsae.flux[as.Date(lh$timeBgn) == days$Day[j]]
      ## Require >50% completeness
      if(length(day) > 24){
        days$LH.max[j] = quantile(day, 0.9, na.rm = TRUE)
      }else{
        days$LH.max[j] = NA
      }
    }
    points(as.POSIXct(days$Day) + 3600 * 12, days$LH.max, pch = 21, bg = "red")
    
    ## LH flux variance for preceding 60 days, require 10 days
    if(sum(!is.na(days$LH.max)) >= 10){
      bouts$lhVar[i] = var(days$LH.max, na.rm = TRUE)
    }else{
      bouts$lhVar[i] = NA
    }
    
    ## Mean max LH flux for preceding 10 days
    days = days[days$Day > as.Date(bouts$Date[i]) - 10,]
    if(nrow(days) > 0){
      if(sum(!is.na(days$LH.max)) > 0){
        bouts$lh10[i] = mean(days$LH.max, na.rm = TRUE)
      }else{
        bouts$lh10[i] = NA
      }
    }else{
      bouts$lh10[i] = NA
    }    
  }else{
    bouts$lhVar[i] = NA
    bouts$lh10[i] = NA
  }
}

write.csv(bouts, "out/Bouts.csv", row.names = FALSE)
