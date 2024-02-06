library(neonUtilities)

p = read.csv("out/plantsNoLow.csv")
bouts = unique(p$Bout)
dts = p$Collection_Date[match(bouts, p$Bout)]
sts = substr(bouts, 1, 4)
bouts = data.frame("Bout" = bouts, "Date" = as.POSIXct(dts), "Site" = sts)
bouts$Precip = rep(0)

for(i in 1:nrow(bouts)){
  pre = loadByProduct("DP1.00006.001", bouts$Site[i], 
                      as.character(bouts$Date[i] - 30 * 3600 * 24), 
                      as.character(bouts$Date[i]), check.size = FALSE)
  if(sum(!(is.na(pre$PRIPRE_30min$priPrecipBulk))) > 0){
    plot(as.POSIXct(pre$PRIPRE_30min$startDateTime), pre$PRIPRE_30min$priPrecipBulk,
         type = "l", main = bouts$Bout[i])
    bouts$Precip[i] = sum(pre$PRIPRE_30min$priPrecipBulk, na.rm = TRUE)
  }else{
    bouts$Precip[i] = NA
  }
}

write.csv(bouts, "out/Bouts.csv")
