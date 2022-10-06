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

r <- raster::getData("worldclim", var="bio", res=2.5)

mat <- r[[1]]
max <- r[[5]]
map <- r[[12]]


meansoil <- raster("raw data/soils/SBIO1_Annual_Mean_Temperature_5_15cm.tif")
plot(meansoil)
rangesoil <- raster("raw data/soils/SBIO7_Temperature_Annual_Range_5_15cm.tif")

#resample soil rasters so they match the climate 

meanresample <- resample(meansoil, mat, method = "bilinear")
rangesoilresample <- resample(rangesoil, mat, method = "bilinear")

envs <- raster::stack(mat, max, map, meanresample, rangesoilresample)
crs <- st_crs(envs)


plot(rangesoilresample)


#crop to north america - this could be better

#left side, right side, bottoms, top
e <- as(extent(-140, -70, 7, 50), 'SpatialPolygons')
#plot(e, add = TRUE)

crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"

envs <- crop(envs, e)
plot(envs[[1]])

## cropped!!

plot(envs[[5]])




#let's make some functions instead of doing this at least 11 times

#cleans data
cleangabi <- function(x, df, env.rast){
  df <- filter(df, valid_species_name == x)
  df <- filter(df, dubious != "Dubious")
  df <- filter(df, dec_lat != "")
  print(unique(df$country))
  df_sf <- st_as_sf(df, coords = c("dec_long", "dec_lat"), crs= st_crs(env.rast))
}

#thin points to one per raster cell
thin_occ <- function(x, env.rast){
  occs.cells <- raster::extract(env.rast[[1]], x, cellnumbers = TRUE)
  occs.cellDups <- duplicated(occs.cells[,1])
  occs <- x
  occs <- occs[!occs.cellDups,]
  occs <- st_transform(occs, crs = st_crs(env.rast), type = "ellps")
}

#generates a convex hull around points, buffers it and then chooses background points from within the hull
generate_bg <- function(x, env.rast){
  hull <- st_convex_hull(st_union(x))
  hull.buf <- sf::st_buffer(hull, dist = 100000)
  hull.buf <- st_sf(hull.buf)
  envs.bg <- raster::mask(env.rast, hull.buf)
  bg <- dismo::randomPoints(envs.bg[[1]], n = 10000) %>% as.data.frame()
}


# run the SDM and extract prediction raster of best model AUC
sdm_seq <- function(occ, env.rast, bgpoints){
  occs <- as.data.frame(st_coordinates(occ))
  colnames(bgpoints) <- colnames(occs)
  e.mx <- ENMevaluate(occs, env.rast, bgpoints, 
                      algorithm = 'maxnet', partitions = 'randomkfold', 
                      tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:5))
  res <- eval.results(e.mx)
  opt.aicc <- res %>% filter(delta.AICc == 0)
  print(opt.aicc)
  opt.seq <- res %>% 
    filter(or.10p.avg == min(or.10p.avg)) %>% 
    filter(auc.val.avg == max(auc.val.avg))
  print(opt.seq)
  preds <- eval.predictions(e.mx)[[opt.aicc$tune.args]]
}

#put the functions together into one!

predraster <- function(x, df, env.rast){
  cleanant <- cleangabi(x, df, env.rast)
  occ <- thin_occ(cleanant, env.rast)
  b <- generate_bg(occ, env.rast)
  sdm <- sdm_seq(occ, env.rast, b)
}


#need to make decisions about partition and how to choose models
#can choose AICc model, and then report AUC in paper



## Solenopsis xyloni


xyloni <- predraster("Solenopsis.xyloni", ants, envs)
plot(xyloni)


## Pheidole hyatti

pheidole <- predraster("Pheidole.hyatti", ants, envs)
plot(pheidole)


dorybi <- predraster("Dorymyrmex.bicolor", ants, envs)

plot(dorybi)

doryin <- predraster("Dorymyrmex.insanus", ants, envs)

plot(doryin)

cypho <- predraster("Cyphomyrmex.wheeleri", ants, envs)
myrmeco <- predraster("Myrmecocystus.kennedyi", ants, envs)
pogo <- predraster("Pogonomyrmex.hoelldobleri", ants, envs)
temno <- predraster("Temnothorax.andrei", ants, envs)
forel <- predraster("Forelius.pruinosus", ants, envs)



messandrei <- predraster("Veromessor.andrei", ants, envs)
messper <- predraster("Veromessor.pergandei", ants, envs)

#Forelius is found in Cuba, others are found in Panama, Nicaragua
#the south american dorymyrmex insanus are likely different species
#Restrict to the extent of N


dorybi@data@names <- "Dorymyrmex bicolor"
doryin@data@names <- "Dorymyrmex insanus"
pheidole@data@names <- "Pheidole hyatti"
xyloni@data@names <- "Solenopsis xyloni"
cypho@data@names <- "Cyphomyrmex wheeleri"
myrmeco@data@names <- "Myrmecocystus"
pogo@data@names <- "Pogonomyrmex californicus"
temno@data@names <- "Temnothorax andrei"
forel@data@names <- "Forelius pruinosis"
messandrei@data@names <- "Messor andrei"
messper@data@names <- "Messor pergandei"

two <- raster::stack(dorybi, doryin, pheidole, xyloni, cypho, myrmeco, pogo, temno, forel, messandrei, messper)
#let's save the prediction raster stack too

saveRDS(two, file = "Clean Data/objects/predictionrasters.rds")

overlap <- calc.niche.overlap(two, overlapStat = "D")
overlap

library(Matrix)
over <- overlap

over2 <- forceSymmetric(over, "L")
over2 <- as.matrix(over2)



#bring in species-level trait data


traits.sp <- traits.sp %>%
  select(.,-sd) %>%
  pivot_wider(names_from = Trait, values_from = mean) %>%
  as.data.frame(traits.sp)

traits.sp <- select(traits.sp, 1, 3, 5, 8, 9, 11, 13, 14)
traits.sp$X.1 <- gsub(" ", ".", traits.sp$X.1)

species.name <- traits.sp$X.1

#gower dissimilarity?
gow <- gowdis(traits.sp)  
gow <- as.matrix(gow)
row.names(gow) <- species.name
euc <- dist(traits.sp, "euclidean")
size <- dist(traits.sp$Femur.w, "euclidean")

over2 <- as.matrix(over2[order(row.names(over2)),])
over2 <- over2[,order(colnames(over2))]

overdis <- 1 - over2

saveRDS(overdis, file = "Clean Data/objects/overlap_dissim.rds")

#save overlap as an object

mantel(gow, overdis)
