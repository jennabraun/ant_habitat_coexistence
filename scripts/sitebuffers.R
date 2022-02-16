#draw a buffer around the studysites and export a shapefile

library(dplyr)
library(sf)

sites <- read.csv("Clean Data/sites_joined.csv")
site_sf <- st_as_sf(sites, coords = c("Long", "Lat"), crs = 4326)
plot(site_sf)
site_buf <- st_buffer(site_sf, dist = 50)
plot(site_buf)
st_write(site_sf, "Clean Data/sites.shp")
