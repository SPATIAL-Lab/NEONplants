library(SIBER)
source("code/multiLikeOverlap.R")

# Using source water medians ----

# CLBJ6
## Prep data
load(paste0("out/CLBJ_tower/CLBJ6.rda"))
ps = data.frame("iso1" = numeric(), "iso2" = numeric(), "group" = character(), 
                "community" = numeric())
for(j in seq_along(smix)){
  psa = data.frame("iso1" = smix[[j]]$summary["mixture_d2H", "50%"], 
                   "iso2" = smix[[j]]$summary["mixture_d18O", "50%"], 
                   "group" = names(smix)[j], 
                   "community" = 1)
  ps = rbind(ps, psa)
}

spf = table(ps$group)
spi = names(spf)[spf > 3]

ps = ps[ps$group %in% spi,]

## SIBER object
ps.t.CLBJ6 = suppressWarnings(createSiberObject(ps))

## Overlap
mlo.CLBJ6 = multiLikOverlap(ps.t.CLBJ6$all.groups, ps.t.CLBJ6, n = 500)

# HARV5
## Prep data
load(paste0("out/HARV_tower/HARV5.rda"))
ps = data.frame("iso1" = numeric(), "iso2" = numeric(), "group" = character(), 
                "community" = numeric())
for(j in seq_along(smix)){
  psa = data.frame("iso1" = smix[[j]]$summary["mixture_d2H", "50%"], 
                   "iso2" = smix[[j]]$summary["mixture_d18O", "50%"], 
                   "group" = names(smix)[j], 
                   "community" = 1)
  ps = rbind(ps, psa)
}

spf = table(ps$group)
spi = names(spf)[spf > 3]

ps = ps[ps$group %in% spi,]

## SIBER object
ps.t.HARV5 = suppressWarnings(createSiberObject(ps))

## Overlap
mlo.HARV5 = multiLikOverlap(ps.t.HARV5$all.groups, ps.t.HARV5, n = 500)

## Plotting args
group.ellipses.args = list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)

# Plot CLBJ6 and HARV5
png("out/mloPlot1.png", width = 9, height = 4, units = "in", res = 600)
layout(matrix(c(1, 2), ncol = 2))
par(mar = c(5,5,3,1))
plotSiberObject(ps.t.CLBJ6,                   
                ax.pad = 5, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = FALSE, group.hull.args,
                bty = "L", x.limits = c(-4, 4),
                iso.order = c(2, 1),
                xlab = expression(delta^18*O~'\u2030'),
                ylab = expression(delta^2*H~'\u2030'))
box()
taxa = lapply(ps.t.CLBJ6$all.groups, function(x) bquote(italic(.(x))))
legend("topleft", as.expression(taxa), cex = 0.8,
       col = seq(length(ps.t.CLBJ6$all.groups)), lty = 1, bty = "n")
title(paste("CLBJ6: Overlap = ", round(mlo.CLBJ6["Fraction Overlap"], 2)))

plotSiberObject(ps.t.HARV5,                   
                ax.pad = 5, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = FALSE, group.hull.args,
                bty = "L", x.limits = c(-10, -6),
                iso.order = c(2, 1),
                xlab = expression(delta^18*O~'\u2030'),
                ylab = expression(delta^2*H~'\u2030'))
box()
taxa = lapply(ps.t.HARV5$all.groups, function(x) bquote(italic(.(x))))
legend("topleft", as.expression(taxa), cex = 0.8,
       col = seq(length(ps.t.HARV5$all.groups)), lty = 1, bty = "n")
title(paste("HARV5: Overlap = ", round(mlo.HARV5["Fraction Overlap"], 2)))

dev.off()

# UNDE7
## Prep data
load(paste0("out/UNDE_tower/UNDE7.rda"))
ps = data.frame("iso1" = numeric(), "iso2" = numeric(), "group" = character(), 
                "community" = numeric())
for(j in seq_along(smix)){
  psa = data.frame("iso1" = smix[[j]]$summary["mixture_d2H", "50%"], 
                   "iso2" = smix[[j]]$summary["mixture_d18O", "50%"], 
                   "group" = names(smix)[j], 
                   "community" = 1)
  ps = rbind(ps, psa)
}

spf = table(ps$group)
spi = names(spf)[spf > 3]

ps = ps[ps$group %in% spi,]

## SIBER object
ps.t.UNDE7 = suppressWarnings(createSiberObject(ps))

## Overlap
mlo.UNDE7 = multiLikOverlap(ps.t.UNDE7$all.groups, ps.t.UNDE7, n = 500)

## Plot
png("out/mloPlot2.png", width = 4.5, height = 4, units = "in", res = 600)
par(mar = c(5,5,3,1))
plotSiberObject(ps.t.UNDE7,                   
                ax.pad = 5, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = FALSE, group.hull.args,
                bty = "L", x.limits = c(-11, -7),
                y.limits = c(-80, -50), iso.order = c(2, 1),
                xlab = expression(delta^18*O~'\u2030'),
                ylab = expression(delta^2*H~'\u2030'))
box()
taxa = lapply(ps.t.UNDE7$all.groups, function(x) bquote(italic(.(x))))
legend("topleft", as.expression(taxa), cex = 0.8,
       col = seq(length(ps.t.UNDE7$all.groups)), lty = 1, bty = "n")
title(paste("UNDE7: Overlap = ", round(mlo.UNDE7["Fraction Overlap"], 2)))

dev.off()