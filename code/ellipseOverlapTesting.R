library(SIBER)

# Prep ####
## Read data
p = read.csv("data/plants.csv")

## Get bits for SIBER
ps = data.frame("iso1" = p$d2H, "iso2" = p$d18O, "group" = p$Species, 
                "community" = p$Bout)

## Make sequential values for community
bouts = unique(ps$community)
ps$community = match(ps$community, bouts)

# Try with one bout for starters ####

## SIBER object
ps.t = createSiberObject(ps[ps$community == 2,])

## Plotting args
group.ellipses.args = list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)

## Plot for fun
par(mar = c(5,5,1,1))
plotSiberObject(ps.t,                   
                ax.pad = 5, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = FALSE, group.hull.args,
                bty = "L",
                iso.order = c(2, 1),
                xlab = expression(delta^18*O~'\u2030'),
                ylab = expression(delta^2*H~'\u2030'))
legend("topleft", ps.t$all.groups, col = seq(length(ps.t$all.groups)), lty = 1)

## Max likelihood overlap
e95 = maxLikOverlap("2.Abies lesiocarpa", "2.Carex rupestris", ps.t, p.interval = 0.95)
e95["overlap"] / (e95["area.1"] + e95["area.2"] - e95["overlap"])

## Try another
e95 = maxLikOverlap("2.Geum rossii", "2.Carex rupestris", ps.t, p.interval = 0.95)
e95["overlap"] / (e95["area.1"] + e95["area.2"] - e95["overlap"])

## Test the Bayesian version
parms = list(n.iter = 2e4, n.burnin = 1e3, n.thin = 5, n.chains = 3)
priors = list(R = diag(2), k = 2, tau.mu = 1e-3)

## Fit Bayesian ellipses
eB = siberMVN(ps.t, parms, priors)
b95 = bayesianOverlap("2.Abies lesiocarpa", "2.Carex rupestris", eB, draws = 100,
                      p.interval = 0.95)
plot(density((b95["overlap"] / (b95["area1"] + b95["area2"] - b95["overlap"]))[,1], 
             from = 0), main = "")

## One more
b95 = bayesianOverlap("2.Geum rossii", "2.Carex rupestris", eB, draws = 500,
                      p.interval = 0.95)
plot(density((b95["overlap"] / (b95["area1"] + b95["area2"] - b95["overlap"]))[,1], 
             from = 0), main = "")

## Run MLE using my multi version...so fast
source("code/multiLikeOverlap.R")
e95 = multiLikOverlap(ps.t$all.groups, ps.t, n = 500)
