#ant SDMs
#tutorial from vignette
# Load packages -- the order here is important because some pkg functions overwrite others.
library(ENMeval)
library(raster)
library(dplyr)
library(sf)
library(tidyr)

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
cleangabi <- function(name, df, env.rast){
  df <- filter(df, valid_species_name == name)
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
sdm_seq <- function(occ, env.rast, bgpoints, name){
  occs <- as.data.frame(st_coordinates(occ))
  colnames(bgpoints) <- colnames(occs)
  e.mx <- ENMevaluate(occs, env.rast, bgpoints, 
                      algorithm = 'maxnet', partitions = 'randomkfold', 
                      tune.args = list(fc = c("L","LQ","LQH","H"), rm = 1:5))
  res <- eval.results(e.mx)
  opt.aicc <- res %>% filter(delta.AICc == 0)
  print(opt.aicc)
  #opt.seq <- res %>% 
    #filter(or.10p.avg == min(or.10p.avg)) %>% 
    #filter(auc.val.avg == max(auc.val.avg))
  #print(opt.seq)
  preds <- eval.predictions(e.mx)[[opt.aicc$tune.args]]
}


# x <- "test.test"
# paste0("Clean Data/objects/", x, "_sdmresults.rds", sep = "")

#put the functions together into one!

predraster <- function(name, df, env.rast){
  cleanant <- cleangabi(name, df, env.rast)
  occ <- thin_occ(cleanant, env.rast)
  b <- generate_bg(occ, env.rast)
  sdm <- sdm_seq(occ, env.rast, b, name)
}


#need to make decisions about partition and how to choose models
#can choose AICc model, and then report AUC in paper

cyph <- cleangabi("Cyphomyrmex.wheeleri", ants, envs)

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
nocypo <- raster::stack(dorybi, doryin, pheidole, xyloni, myrmeco, pogo, temno, forel, messandrei, messper)
#let's save the prediction raster stack too



saveRDS(two, file = "Clean Data/objects/predictionrasters.rds")
two <- readRDS(file = "Clean Data/objects/predictionrasters.rds")
overlap <- calc.niche.overlap(two, overlapStat = "D")
testover <- 1 - overlap

library(Matrix)
over <- overlap

over2 <- forceSymmetric(over, "L")
over2 <- as.matrix(over2)


#bring in species-level trait data


traits.sp1 <- traits.sp %>%
  dplyr::select(.,-sd) %>%
  pivot_wider(names_from = Trait, values_from = mean) %>%
  as.data.frame(traits.sp)

traits.sp1 <- dplyr::select(traits.sp1, 1, 3, 5, 8, 9, 11, 13, 14)
traits.sp1$X.1 <- gsub(" ", ".", traits.sp1$X.1)

species.name <- traits.sp1$X.1
traits.sp1 <- dplyr::select(traits.sp1, -X.1)
#gower dissimilarity?
gow <- gowdis(traits.sp1)  
gow <- as.matrix(gow)
row.names(gow) <- species.name


over2 <- as.matrix(over2[order(row.names(over2)),])
over2 <- over2[,order(colnames(over2))]

overdis <- 1 - over2

#save overlap as an object
saveRDS(overdis, file = "Clean Data/objects/overlap_dissim.rds")

overdis <- readRDS(file = "Clean Data/objects/overlap_dissim.rds")


gow1 <- as.matrix(gow) 
gow1 <- gow1[lower.tri(gow1)]

overdis1 <- overdis[lower.tri(overdis)]
overdis1 <- as.dist(overdis)
gow1 <- as.matrix(gow1)
mantel(gow, overdis1)


mantel(gow, overdis)
#try again without cyphomyrmex

mdata <- cbind(gow1, overdis1)
mdata <- as.data.frame(mdata)


ggplot(mdata, aes(overdis1, V1)) +
  geom_point() + 
  geom_smooth(method = "lm", colour = "black") + xlab("Climatic Niche Dissimilarity (Complement of Shoener's D)") + 
  ylab("Trait Dissimilarity (Gower Distance)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 13)) +
  geom_label(aes(x = 0.70, y = 0.58),vjust=1, hjust = 0, 
             label = "Mantel r = 0.378, \np = 0.028")





overlapnocyph <- calc.niche.overlap(nocypo, overlapStat = "D")
overlapnocyph <- forceSymmetric(overlapnocyph, "L")
overlapnocyph <- as.matrix(overlapnocyph)


traits.sp.n <- traits.sp %>%
  select(.,-sd) %>%
  filter(X.1 != "Cyphomyrmex wheeleri")%>%
  pivot_wider(names_from = Trait, values_from = mean) %>%
  as.data.frame(traits.sp.n)


traits.sp.n <- select(traits.sp.n, 1, 3, 5, 8, 9, 11, 13, 14)
traits.sp.n$X.1 <- gsub(" ", ".", traits.sp.n$X.1)


gow2 <- gowdis(traits.sp.n)  
gow2 <- as.matrix(gow2)
row.names(gow2) <- species.name
euc2 <- dist(traits.sp.n, "euclidean")
#size <- dist(traits.sp$Femur.w, "euclidean")

overlapnocyph <- as.matrix(overlapnocyph[order(row.names(overlapnocyph)),])
overlapnocyph <- overlapnocyph[,order(colnames(overlapnocyph))]

overdis2 <- 1 - overlapnocyph



mantel(gow2, overdis2)



#extracting model info from rds files
cw <- readRDS("Clean Data/objects/Cyphomyrmex.wheeleri_sdmresults.rds")
res <- eval.results(cw)
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc
eval.variable.importance(cw)[[opt.aicc$tune.args]]
m <- cw@models$fc.H_rm.2
eval.variable.importance(m)
eval.varimp(cw)
modcw <- mod.seq <- eval.models(cw)[[opt.aicc$tune.args]]
modcw$betas
plot(modcw, type = "cloglog")
eval.variable.importance(modcw)
var <- cw@variable.importance

aic.mod <- cw@models[[which(cw@results$delta.AICc==0)]]
eval.variable.importance(aic.mod)
