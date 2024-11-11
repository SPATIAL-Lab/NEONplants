library(terra)
library(tmap)
library(sf)
library(assignR)

sites = vect("out/Sites.shp")

ex = st_bbox(c(xmin = -125, ymin = 25, xmax = -65, ymax = 50), crs = "WGS84")

set_defaults(map_service = "esri", map_type = "world_terrain_base")
b = basemap_terra(ex)

png("out/siteMap.png", 8, 8 * dim(b)[1] / dim(b)[2], units = "in", res = 600)
plot(b, box = TRUE)
plot(project(states, b), add = TRUE, border = "darkgrey")
text(project(sites, b), sites$Site, halo = TRUE)
dev.off()
