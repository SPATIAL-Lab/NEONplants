library(isoWater)
library(assignR)
library(sp)
library(raster)

#Get soil and plant data from wiDB
d = wiDB_data(projects = '00384')
d = d$data

#Subset data to get one site
unique(d$Site_ID)
d.wood = d[d$Site_ID == "WOOD_tower",]
