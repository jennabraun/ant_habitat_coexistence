Distribution <- as.matrix(rbind(c(1,1,0,0),c(1,0,0,0),c(0,0,1,1),c(0,0,1,0)))
T <- c(1,2,5,8)
E <- c(10,12,100,112)
n.species <- ncol(Distribution)
n.sites <- nrow(Distribution)
Dist.stacked <- as.vector(Distribution)
E.stacked <- rep(1, n.species) %x% scale(E)
T.stacked <- scale(T) %x% rep(1, n.sites)
predictor <- E.stacked * T.stacked
model <- glm(Dist.stacked ~ predictor,family=binomial(link=logit))
summary(model)
predictor2 <- T %x% E

library(ade4)
E <- aravo$env[,c("Aspect","Slope","Snow")]
E <- as.matrix(E)
T <- as.matrix(aravo$trait)
TE <- kronecker(T, E, make.dimnames = TRUE)
dim(TE)
colnames(TE)



n.species <- ncol(Distribution)
n.sites <- nrow(Distribution)
Dist.stacked <- as.vector(Distribution)
E.stacked <- rep(1, n.species) %x% E
T.stacked <- T %x% rep(1, n.sites)  

View(cbind(Dist.stacked,E.stacked,T.stacked))


Dist.stacked <- as.vector(as.matrix(aravo$spe))
ColNames <- colnames(TE)
TE <- scale(TE) # so that slopes can be compared directly to one another
colnames(TE) <- ColNames
model.bilinear.poisson <- glm(Dist.stacked ~ TE,family="poisson")
library(MASS)
model.bilinear.negBinom <- glm.nb(Dist.stacked  ~ TE)


set.seed(125)
E <- rnorm(n.communities)
T <- rnorm(n.species)
Dist <- generate_communities(tolerance = 1.5, E, T, n.species, n.communities)
Distribution <- Dist$Dist.matrix
T.false <- rnorm(n.species)
TraitEnv.res <- TraitEnvCor(Distribution,E,T.false, Chessel = FALSE)
TraitEnv.res["Fourthcorner"]

Distribution <- as.matrix(rbind(c(1,1,0,0),c(1,0,0,0),c(0,0,1,1),c(0,0,1,0)))
Distribution
T <- c(1,2,5,8)
E <- c(10,12,100,112)
source("UtilityFunctions.R")
TraitEnvCor(Distribution,E,T)["Fourthcorner"]



E <- aravo$env[,c("Aspect","Slope","Snow")]
E <- as.matrix(E)
T <- as.matrix(aravo$trait)       
TE <- kronecker(T, E, make.dimnames = TRUE)
dim(TE)
colnames(TE)


Dist.stacked <- as.vector(as.matrix(aravo$spe))
ColNames <- colnames(TE)
Dist.stacked
TE <- scale(TE) # so that slopes can be compared directly to one another
colnames(TE) <- ColNames
model.bilinear.poisson <- glm(Dist.stacked ~ TE,family="poisson")
library(MASS)
model.bilinear.negBinom <- glm.nb(Dist.stacked  ~ TE)
