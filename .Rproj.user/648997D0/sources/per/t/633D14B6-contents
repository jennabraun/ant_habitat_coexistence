library(dplyr)
library(performance)
library(ggplot2)

#quick PCA to display differences in vegetation
veg <- read.csv("raw data/vegetation.csv")

veg$sampleID <- paste(veg$Site, veg$month, veg$Microsite, veg$Rep)
veg <- rename(veg, "Dry Veg" = dry.veg.percent, "Green Veg" = green.veg.percent, "Woody Debris" = woody.debris.percent, "Rocks" = rocks.percent, "Bare Ground" = bare.ground.percent)
veg$shrub.x[is.na(veg$shrub.x)] <- 0
veg$shrub.y[is.na(veg$shrub.y)] <- 0
veg$shrub.z[is.na(veg$shrub.z)] <- 0
veg <- na.omit(veg)
veg <- select(veg, -16,-17, -18)
var <- select(veg, sampleID, 9:13, 15)
cov <- select(veg, sampleID, 1, 2, 4)
row.names(var) <- var$sampleID
row.names(cov) <- cov$sampleID
all.equal(rownames(var), rownames(cov))

#drop sample id from var
var <- select(var, -sampleID)




ggplot(veg, aes(Microsite, dry.veg.percent)) + geom_violin() + geom_boxplot(width = 0.05) 
  





pca1<- prcomp(var)
summary(pca1)
pca1

library(ggfortify)
autoplot(pca1, data = cov, colour = "Microsite", loadings = TRUE, loadings.label = TRUE)



m1 <- glm(dry.veg.percent ~ Microsite, data = veg, family = poisson)
summary(m1)
check_overdispersion(m1)
m2 <- MASS::glm.nb(dry.veg.percent ~ Microsite, data = veg)
summary(m2)
