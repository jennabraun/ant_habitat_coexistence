library(raster)
library(sp)
library(sf)
library(dplyr)

sites <- read.csv("Clean data/studysites.csv")
site <- sites$site
coords <- data.frame(x=sites$long,y=sites$lat)

points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

#northern sites
r <- getData("worldclim",var="bio",res=0.5, lon = sites[1,3], lat = sites[1,2])
r <- r[[c(1,12,5)]]
names(r) <- c("Temp","Prec", "Max")
values <- extract(r,points)
north <- cbind.data.frame(coordinates(points),values)
north
north$site <- site

#southern sites
r <- getData("worldclim",var="bio",res=0.5, lon = sites[5,3], lat = sites[5,2])
r <- r[[c(1,12,5)]]
names(r) <- c("Temp","Prec", "Max")
values <- extract(r,points)
south <- cbind.data.frame(coordinates(points),values)
south
south$site <- site




#put together
north <- filter(north, Temp > 0)
south <- filter(south, Temp > 0)
sites <- rbind(north, south)

sites <- mutate(sites, Temp = Temp/10)
sites <- mutate(sites, Max = Max/10)
sites <- mutate(sites, arid = Prec/(Temp+10))
write.csv(sites, "Clean data/sites_worldclim.csv")
