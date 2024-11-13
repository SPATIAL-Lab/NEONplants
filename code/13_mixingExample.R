library(MCMCpack)

# Load data
## CLBJ6 mixing results
load("out/CLBJ_tower/CLBJ6.rda")
load("out/sourceStats.rda")

## Pick a xylem sample
i = 7
p = read.csv("out/plantsNoLow.csv")
p = p[p$Bout == "CLBJ6",]
p = p[i, ]
smix = smix[[i]]
KLd = round(source.stats$CLBJ6$KLd[i], 1)

## Prior
prior = rdirichlet(1e6, rep(0.3, 5))
prior = density(prior[, 1])

png("out/mixingExample.png", width = 6, height = 4, units = "in", res = 600)
ytop = -2
par(mar = c(5, 5, 3, 1))
plot(sources$O, sources$H, ylim = c(-36, ytop), xlim = c(-5, 3), axes = FALSE,
     xlab = expression(delta^{18}*"O"), ylab = expression(delta^2*"H"))
points(smix$results$mixture_d18O, smix$results$mixture_d2H,
       pch = ".")
points(sources$O, sources$H, cex = 2, lwd = 2, pch = 21, bg = "salmon2")
points(p$d18O.oc, p$d2H, pch = 21, cex = 2, lwd = 2, bg = "seagreen")
for(i in 1:nrow(sources)){
  d = density(smix$results[, 2 + i], from = 0, to = 1)
  polygon(c(sources$O[i] + d$x - 0.5, sources$O[i] + 0.5, sources$O[i] - 0.5), 
          c(ytop + d$y, ytop, ytop), col = "salmon2", xpd = NA)
  lines(sources$O[i] + d$x - 0.5, ytop + d$y, xpd = NA, lwd = 2)
  lines(sources$O[i] + prior$x - 0.5, ytop + prior$y, xpd = NA, lty = 3)
}
lines(c(min(sources$O) - 0.5, max(sources$O) + 0.5), rep(ytop, 2))
text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]),
     par("usr")[4] - 0.1 * diff(par("usr")[3:4]),
     paste("Kullback-Leibler divergence =", KLd, "bit"), pos = 4)
axis(1)
axis(2)

dev.off()
