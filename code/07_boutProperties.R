library(neonUtilities)

# Prep ----
p = read.csv("out/plantsNoLow.csv")
bouts = unique(p$Bout)
dts = p$Collection_Date[match(bouts, p$Bout)]
sts = substr(bouts, 1, 4)
bouts = data.frame("Bout" = bouts, "Date" = as.POSIXct(dts), "Site" = sts)

# Precipitation ----
bouts$Precip = rep(0)

## Cycle through bouts
for(i in 1:nrow(bouts)){
  pre = loadByProduct("DP1.00006.001", bouts$Site[i], 
                      as.character(bouts$Date[i] - 30 * 3600 * 24), 
                      as.character(bouts$Date[i]), check.size = FALSE)
  pre = pre$PRIPRE_30min
  pre = pre[pre$startDateTime > bouts$Date[i] - 30 * 3600 * 24,]
  pre = pre[pre$startDateTime < bouts$Date[i],]
  if(sum(!(is.na(pre$priPrecipBulk))) > 0){
    plot(as.POSIXct(pre$startDateTime), pre$priPrecipBulk,
         type = "l", main = bouts$Bout[i])
    bouts$Precip[i] = sum(pre$priPrecipBulk, na.rm = TRUE)
  }else{
    bouts$Precip[i] = NA
  }
}

# Soil moisture ----
bouts$VWC1 = bouts$VWC3 = rep(0)

## Cycle through bouts
for(i in 1:nrow(bouts)){
  vwc = loadByProduct("DP1.00094.001", bouts$Site[i], 
                      as.character(bouts$Date[i] - 30 * 3600 * 24), 
                      as.character(bouts$Date[i]), check.size = FALSE)
  vwc = vwc$SWS_30_minute
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
    
    bouts$VWC1[i] = mean(vwc1$VSWCMean, na.rm = TRUE)
    bouts$VWC3[i] = mean(vwc3$VSWCMean, na.rm = TRUE)
    bouts$VWC5[i] = mean(vwc3$VSWCMean, na.rm = TRUE)
  }else{
    bouts$VWC1[i] = bouts$VWC3[i] = bouts$VWC5[i] = NA
  }
}

write.csv(bouts, "out/Bouts.csv")
