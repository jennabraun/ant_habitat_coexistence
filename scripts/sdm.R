#ant SDMs
#tutorial from vignette
# Load packages -- the order here is important because some pkg functions overwrite others.
library(ENMeval)
library(raster)
library(dplyr)
library(sf)
# Set a random seed in order to be able to reproduce this analysis.
set.seed(48)


#read in GABI ants data
ants <- read.csv("raw data/GABI_Data_Release1.0_18012020.csv")
fire <- filter(ants, valid_species_name == "Solenopsis.xyloni")
fire <- filter(fire, dubious != "Dubious")
#I'm removing the Bahamas population
fire <- filter(fire, country != "Bahamas")
unique(fire$country)


write.csv(fire, "fireants.csv")


fire <- filter(fire, dec_lat != "")
crs <- st_crs(envs)
fire_sf <- st_as_sf(fire, coords = c("dec_long", "dec_lat"), crs= st_crs(sh))

st_write(fire_sf, "fire.shp")
occs.cells <- raster::extract(envs[[1]], fire_sf, cellnumbers = TRUE)
occs.cellDups <- duplicated(occs.cells[,1])
occs <- fire_sf
occs <- occs[!occs.cellDups,]
occs <- st_transform(occs, crs = st_crs(envs), type = "ellps")


plot(envs[[1]], main="Mean annual temperature")
plot(occs, add = TRUE)
#can compare outputs with GBIF too

?sf
sh <- st_read("raw data/Bentity2_shapefile_fullres.shp")
plot(sh, add = TRUE)


st_crs(sh)
st_crs(envs)
st_crs(occs)
# You can search online databases like GBIF using the spocc package (commented below),
# but here we will load in some pre-downloaded data.
# bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=300, has_coords=TRUE)
# occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
occs <- readRDS("bvariegatus.rds")

# Removing occurrences that have the same coordinates is good practice to
# avoid pseudoreplication.
occs <- occs[!duplicated(occs),]


envs.files <- list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), 
                         pattern='grd', full.names=TRUE)
# Read the raster files into a RasterStack.
# These variables represent 8 bioclimatic variables and one categorical variable "biome".
# Find the descriptions of the bioclimatic variables here: 
# https://www.worldclim.org/data/bioclim.html
envs <- raster::stack(envs.files)
# The biome raster has some NAs for cells that have values in the other rasters.
# Let's mask all rasters to biome to change the value of these cells to NA for all rasters.
# ENMeval will do this automatically, but let's do it here to avoid the warning message later.
# We change back from a RasterBrick to RasterStack because of issues with assigning 
# factor rasters for RasterBricks.
envs <- raster::mask(envs, envs[[9]]) %>% raster::stack()
# Make sure to declare the categorical variable as a factor
envs$biome <- raster::as.factor(envs$biome)
# Let's now remove occurrences that are cell duplicates -- these are
# occurrences that share a grid cell in the predictor variable rasters.
# Although Maxent does this by default, keep in mind that for other algorithms you may
# or may not want to do this based on the aims of your study.
# Another way to space occurrence records a defined distance from each other to avoid
# spatial autocorrelation is with spatial thinning (Aiello-Lammens et al. 2015).
occs.cells <- raster::extract(envs[[1]], occs, cellnumbers = TRUE)

occs.cellDups <- duplicated(occs.cells[,1])
occs <- occs[!occs.cellDups,]

# Plot first raster in the stack, the mean annual temperature.
plot(envs[[1]], main="Mean annual temperature")

# Add points for all the occurrence points onto the raster.
points(occs)

# There are some points east of the Amazon River.
# Suppose we know this is a population that we don't want to include in the model.
# We can remove these points from the analysis by subsetting the occurrences by 
# latitude and longitude.
occs <- filter(occs, latitude > -20, longitude < -45)

# Plot the subsetted occurrences to make sure we filtered correctly.
points(occs, col = 'red')