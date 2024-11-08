pr = read.csv("out/plantWDex.csv")
p = read.csv("out/plantsNoLow.csv")
s = read.csv("data/soils.csv")

# All, unscreened
xlim = range(c(pr$d18O, s$d18O))
ylim = range(c(pr$d2H, s$d2H))

png("out/allData.png", 6, 6, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
plot(pr$d18O, pr$d2H, xlim = xlim, ylim = ylim, pch = 20, cex = 0.5,
     xlab = expression(delta^{18}*"O (VSMOW)"),
     ylab = expression(delta^2*"H (VSMOW)"))
abline(10, 8, lwd = 2)
points(s$d18O, s$d2H, pch = 21, bg = "salmon2")
points(pr$d18O, pr$d2H, pch = 21, bg = "seagreen")
dev.off()

# All, screened
xlim = range(c(p$d18O.oc, s$d18O))
ylim = range(c(p$d2H, s$d2H))

png("out/allDataScreened.png", 6, 6, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
plot(p$d18O.oc, p$d2H, xlim = xlim, ylim = ylim, pch = 20, cex = 0.5,
     xlab = expression(delta^{18}*"O (VSMOW)"),
     ylab = expression(delta^2*"H (VSMOW)"))
abline(10, 8, lwd = 2)
points(s$d18O, s$d2H, pch = 21, bg = "salmon2")
points(p$d18O.oc, p$d2H, pch = 21, bg = "seagreen")
dev.off()
