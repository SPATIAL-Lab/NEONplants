# overlap for n group ellipses 
library(terra)

multiLikOverlap = function(groups, siber.object, 
                          p.interval = 0.95, n = 100,
                          do.plot = FALSE) {
  
  # ----------------------------------------------------------------------------

  coords = list()

  for(i in seq_along(groups)){
    coords[[i]] = addEllipse(siber.object$ML.mu[[1]][ , , groups[i]],
                           siber.object$ML.cov[[1]][ , , groups[i]],
                           m = siber.object$sample.sizes[1, groups[[i]]],
                           small.sample = TRUE,
                           n = n,
                           p.interval = p.interval,
                           ci.mean = FALSE,
                           do.plot = FALSE)

  }
  
  polys = vect(coords, type = "poly")
  polys = rbind(polys, aggregate(polys))
  areas = suppressWarnings(expanse(polys))
  
  areas = append(areas, 1 - (areas[length(areas)] / sum(areas[-length(areas)])))
  
  names(areas) = c(groups, "Combined", "Fraction Overlap")

  return(areas)
}