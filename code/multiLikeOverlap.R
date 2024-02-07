# overlap for n group ellipses 
multiLikOverlap = function(groups, siber.object, 
                          p.interval = 0.95, n = 100) {
  
  library(terra)
  
  if(length(groups) < 2){
    stop("Need at least 2 groups")
  }
  
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
  areas = suppressWarnings(expanse(polys))
  
  sumarea = sumoarea = 0
  for(i in 1:(length(groups) - 1)){
    for(j in (i + 1):length(groups)){
      sumarea = sumarea + sum(suppressWarnings(expanse(polys[c(i, j)])))
      oarea = suppressWarnings(expanse(intersect(polys[i], polys[j])))
      if(length(oarea) == 0){oarea = 0}
      sumoarea = sumoarea + oarea
    }
  }
  
  areas = append(areas, sumoarea / sumarea)
  
  names(areas) = c(groups, "Fraction Overlap")

  return(areas)
}