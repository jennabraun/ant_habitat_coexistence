#exploring the TPD packages


library(TPD)
head(iris)
traits_iris <- iris[, c("Sepal.Length", "Sepal.Width")]
sp_iris <- iris$Species
TPDs_iris <- TPDs(species = sp_iris, traits_iris)

head(TPDs_iris$data$evaluation_grid)
nrow(TPDs_iris$data$evaluation_grid)
names(TPDs_iris$TPDs)
sum(TPDs_iris$TPDs$setosa)
pops_iris <- rep(c(rep("Site.1", 25), rep("Site.2", 25)), 3)
TPDs_pops_iris <- TPDs (species = sp_iris, traits = traits_iris, samples = pops_iris)




library(FD)
dummy


library(picante)
rep.is <- replicate(5, randomizeMatrix(wide.pop, null.model = "independentswap"))

numberReps <- 100
#Lets create a matrix to store results from each iteration (one column per iteration)
resultsRandom <- matrix(NA, nrow = nrow(wide.pop), ncol = numberReps,
                        dimnames = list(rownames(wide.pop),
                                        paste0("Sim.", 1:numberReps)))
for(rep in 1:numberReps){
  randomizedFDis <- randomizeMatrix(samp = wide.pop, null.model = "independentswap")
  simFDis <- dbFD(trait.pop, randomizedFDis, w.abun = TRUE)$FDis
  resultsRandom[, rep] <- simFDis
}

obsFDis <- dbFD(trait.pop, wide.pop, w.abun = TRUE)$FDis
meanNull3 <- rowMeans(resultsRandom)
ES3 <- obsFDis - meanNull3
sdNull3 <- apply(resultsRandom, 1, sd)
SES3 <- ES3 / sdNull3
SES_dis <- data.frame(ES3)
