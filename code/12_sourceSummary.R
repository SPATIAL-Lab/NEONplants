library(MCMCpack)

# Get bouts and sites
p = read.csv("out/plantsNoLow.csv")

# Unique sites
bouts = unique(p$Bout)
sites = p$Site_ID[match(bouts, p$Bout)]

## Index for source depths
depths = c("0-10 cm", "10-20 cm", "20-40 cm", "40+ cm", "Groundwater")

## Function to summarize source depths
sprod = function(fracs, dind){
  dind = dind[!is.na(dind)]
  wfracs = fracs * dind
  return(sum(wfracs))
}

## Prior distributions on source fractions
priors = list()
bw = 0.05
alpha = 0.3
for(i in 1:5){
  priors[[i]] = density(rdirichlet(1e6, rep(alpha, i))[, 1], from = 0, to = 1, 
                        bw = bw)
  priors[[i]]$y = priors[[i]]$y / sum(priors[[i]]$y)
}

## Space for results
source.stats = list()

## Loop through sites
for(i in seq_along(bouts)){
  ## Space for site
  source.stats[[i]] = list()
  names(source.stats)[i] = bouts[i]

  ## Load smix and sources
  load(file.path("out", sites[i], paste0(bouts[i], ".rda")))

  ## Calculate and save source ranges
  srange = c(diff(range(sources$H)), diff(range(sources$O)))
  source.stats[[i]][[1]] = srange

  ## Space for results
  sdepths = matrix(nrow = length(smix), ncol = 7500)
  taxa = character()
  kld = numeric()
  
  ## Loop through Individuals
  for(j in 1:length(smix)){
    dind = match(names(smix[[j]]$results), depths)
    sdepths[j,] = apply(smix[[j]]$results[, !is.na(dind)], 1, sprod, dind)
    taxa = c(taxa, names(smix)[j])
    
    ## Kullback-Leibler divergence
    ns = sum(!(is.na(dind)))
    if(ns > 1){
      kld = c(kld, 0)
      for(k in 1:(ns - 1)){
        post = density(smix[[j]]$results[, 2 + k], from = 0, to = 1, 
                       bw = bw)
        post$y = post$y / sum(post$y)
        
        klds = sum(post$y * log(post$y / priors[[ns]]$y, base = 2), 
                   na.rm = TRUE)
        kld[length(kld)] = kld[length(kld)] + klds
      }
    }
  }
  
  ssd = apply(sdepths, 2, sd)
  
  source.stats[[i]][[2]] = taxa
  source.stats[[i]][[3]] = sdepths
  source.stats[[i]][[4]] = ssd
  source.stats[[i]][[5]] = kld
  
  names(source.stats[[i]]) = c("source.range", "taxa", 
                                    "source.depths", "sd.sd", "KLd")
}

save(source.stats, file = "out/sourceStats.rda")

# Summarize variation in rooting depth estimates for individuals
## Space for results
nr = length(bouts)
ssd.stats = data.frame("Bout" = character(nr), "n" = numeric(nr), 
                       "srangeO" = numeric(nr), "SSD025" = numeric(nr), 
                       "SSD25" = numeric(nr), "SSD50" = numeric(nr), 
                       "SSD75" = numeric(nr), "SSD975" = numeric(nr))

## Loop through bouts
for(i in seq_along(bouts)){
  ssd.stats$Bout[i] = bouts[i]
  ssd.stats$n[i] = length(source.stats[[i]]$taxa)
  ssd.stats$srangeO[i] = source.stats[[i]]$source.range[2]
  if(ssd.stats$n[i] > 1){
    ssd.stats[i, 4:8] = quantile(source.stats[[i]]$sd.sd, 
                                 c(0.025, 0.25, 0.5, 0.75, 0.975))
  } else{
    ssd.stats[i, 4:8] = rep(NA)
  }
}

## Is this dependent on source d18O range?
plot(ssd.stats$srangeO, ssd.stats$SSD50, ylim = c(0.3, 1.3), pch = 21, 
     bg = "seagreen", xlab = expression("Source water "*delta^18*"O range"),
     ylab = c("SD of median source water depth posterior"))
for(i in 1:nrow(ssd.stats)){
  lines(rep(ssd.stats$srangeO[i], 2), c(ssd.stats$SSD25[i], 
                                        ssd.stats$SSD75[i]))
}
summary(lm(SSD50 ~ srangeO, ssd.stats))

## Is this dependent on number of individuals?
plot(ssd.stats$n, ssd.stats$SSD50, ylim = c(0.3, 1.3))
for(i in 1:nrow(ssd.stats)){
  lines(rep(ssd.stats$n[i], 2), c(ssd.stats$SSD25[i], 
                                        ssd.stats$SSD75[i]))
}

write.csv(ssd.stats, "out/ssd.csv", row.names = FALSE)
