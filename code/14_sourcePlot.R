source("code/violin.R")

# Set up
load("out/sourceStats.rda")

# Example plot CLBJ6
ss = source.stats[[match("CLBJ6", names(source.stats))]]
CLBJ6 = violin.prep(ss)

# HARV5
ss = source.stats[[match("HARV5", names(source.stats))]]
HARV5 = violin.prep(ss)

# Plot it
sm = 0.6
png("out/sourceDepths.png", width = 13 * sm, height = 9 * sm,
    units = "in", res = 600)
layout(matrix(c(1, 2), ncol = 2))
violin(CLBJ6, col = "seagreen", lwd = 2, main = "CLBJ6")
violin(HARV5, col = "seagreen", lwd = 2, main = "HARV5")
dev.off()
