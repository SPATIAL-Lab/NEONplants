# Read data
p = read.csv("out/plantsNoLow.csv")
s = read.csv("data/soils.csv")

# Groundwater
gw = read.csv("out/gw.csv")

# CLBJ6, Sept 29 2021
p.CLBJ6 = p[p$Bout == "CLBJ6",]
s.CLBJ6 = s[s$Bout == "CLBJ6",]
gw.CLBJ = gw[gw$Site == "CLBJ_tower",]

# HARV5, July 5, 2021
p.HARV5 = p[p$Bout == "HARV5",]
s.HARV5 = s[s$Bout == "HARV5",]
gw.HARV = gw[gw$Site == "HARV_tower",]

# Plot
png("out/depthPlot.png", width = 9, height = 5, units = "in", res = 600)
layout(matrix(c(1, 2), ncol = 2))

## CLBJ6
par(mar = c(5, 5, 3, 1))
plot(s.CLBJ6$d18O, s.CLBJ6$Depth_meters, ylim = c(0.85, -0.1), axes = FALSE,
     xlim = range(c(s.CLBJ6$d18O, p.CLBJ6$d18O.oc, 
                    gw.CLBJ$d18O.m - gw.CLBJ$d18O.sd)),
     xlab = expression(delta^{18}*"O"), ylab = "Depth (m)", pch = 21,
     cex = 2, lwd = 2, bg = "salmon2", main = "CLBJ6: Sept. 29, 2021")
arrows(gw.CLBJ$d18O.m - gw.CLBJ$d18O.sd, 0.8,
       gw.CLBJ$d18O.m + gw.CLBJ$d18O.sd, 0.8, lwd = 2, length = 0)
points(gw.CLBJ$d18O.m, 0.8, pch = 21, cex = 2, lwd = 2, bg = "salmon2")
axis(1)
axis(2, c(0, 0.3, 0.6))
axis(2, 0.8, "GW")
lines(par("usr")[1:2], c(0, 0), lty = 3)
points(p.CLBJ6$d18O.oc, rep(-0.05, nrow(p.CLBJ6)), pch = 21, cex = 2,
       lwd = 2, bg = "seagreen")

## HARV5
plot(s.HARV5$d18O, s.HARV5$Depth_meters, ylim = c(0.85, -0.1), axes = FALSE,
     xlim = range(c(s.HARV5$d18O, p.HARV5$d18O.oc, 
                    gw.HARV$d18O.m - gw.HARV$d18O.sd)),
     xlab = expression(delta^{18}*"O"), ylab = "Depth (m)", pch = 21,
     cex = 2, lwd = 2, bg = "salmon2", main = "HARV5: July 5, 2021")
arrows(gw.HARV$d18O.m - gw.HARV$d18O.sd, 0.8,
       gw.HARV$d18O.m + gw.HARV$d18O.sd, 0.8, lwd = 2, length = 0)
points(gw.HARV$d18O.m, 0.8, pch = 21, cex = 2, lwd = 2, bg = "salmon2")
axis(1)
axis(2, c(0, 0.3, 0.6))
axis(2, 0.8, "GW")
lines(par("usr")[1:2], c(0, 0), lty = 3)
points(p.HARV5$d18O.oc, rep(-0.05, nrow(p.HARV5)), pch = 21, cex = 2,
       lwd = 2, bg = "seagreen")

dev.off()
