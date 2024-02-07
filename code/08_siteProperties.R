library(prism)
library(terra)
library(httr)
library(jsonlite)

# PRISM ----
prism_set_dl_dir(tempdir())

get_prism_normals("ppt", "4km", annual = TRUE, keepZip = FALSE)
get_prism_normals("vpdmax", "4km", annual = TRUE, keepZip = FALSE)

ppt = rast(file.path(tempdir(), "PRISM_ppt_30yr_normal_4kmM4_annual_bil/PRISM_ppt_30yr_normal_4kmM4_annual_bil.bil"))
vpd = rast(file.path(tempdir(), "PRISM_vpdmax_30yr_normal_4kmM5_annual_bil/PRISM_vpdmax_30yr_normal_4kmM5_annual_bil.bil"))
plot(vpd)

# NEON site data ----
bouts = read.csv("out/Bouts.csv")
sites = data.frame("Site" = unique(bouts$Site), "Lon" = rep(0), "Lat" = rep(0))

for(i in seq_along(sites$Site)){
  loc = GET(paste0("http://data.neonscience.org/api/v0/locations/", sites$Site[i]))  
  loc = content(loc, as = "text")
  loc = fromJSON(loc, simplifyDataFrame = TRUE, flatten = TRUE)
  sites$Lon[i] = loc$data$locationDecimalLongitude
  sites$Lat[i] = loc$data$locationDecimalLatitude
}

sites = vect(sites, geom = c("Lon", "Lat"), crs = "WGS84")
points(sites)

# Extract at sites ----
clim = c(ppt, vpd)
names(clim) = c("MAP", "VPDmax")
sites = cbind(sites, extract(clim, sites, ID = FALSE))
writeVector(sites, "out/Sites.shp")
