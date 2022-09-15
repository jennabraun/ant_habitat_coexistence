# code scraps





#calculate CWM trait values

# This is for abundance
#it's proportional abundances 

site.j <- select(cov, uniID, Site)
wide <- right_join(site.j, wide, by = "uniID")

wide <- select(wide, -X)
#remove the two other solenopsis species with 1 individual
wide <- select(wide, -12, -13)
#wide <- rename(wide, Site = Site.x) %>%
#select(-Site.y)

long <- pivot_longer(wide, 3:13, names_to = "species", values_to = "count")
#get total abundance per trap
long <- group_by(long, uniID) %>% 
  mutate(total = sum(count)) 

long <- long %>% mutate(rel.abun = count/total)

traits.ag <- rename(traits.ag, species = X.1)
traits.ag$species <- gsub(" ", "", traits.ag$species)

#this adds NA traits to 0 abundances - species isn't found at pitfall or site
join.right <- right_join(traits.ag, long, by = c("Site", "species"))

#join.left <- left_join(traits.ag, long, by = c("Site", "species"))

CWM <- mutate(join.right, cwm.prod = mean*rel.abun)

CWM <- group_by(CWM, Site, uniID, Trait) %>% summarise(cwm = sum(cwm.prod))

#can log transform the abundances to look at dominants vs subordinates



#join environmental data to CWM trait dataframe

CWM <- left_join(CWM, cov, by = "uniID")




## NMDS

```{r}

#get rid of zero rows in species data and also in environmental data
zero <- filter(spe.bin,rowSums(spe.bin)!=0)
zeroID <- row.names(zero)
zeroenv <- cov[cov$uniID %in% zeroID, ]
envz <- select(zeroenv,2,  3:9, 12, 15, 17, 23, 29, 34)
#wow so many zeros?

#nm <- metaMDS(zero)

#envfit(nm ~ esi +mean.cover + dry.veg.percent + shrub.site + Microsite + var.cover + veg..height + mean.height + Temp,envz, strata = envz$Site, perm = 999)


# m1 <- cca(zero ~ esi +mean.cover + dry.veg.percent + shrub.site + Microsite + var.cover + veg..height + mean.height + Temp + Site,envz, perm = 999)
# summary(m1)
# anova.cca(m1, by = "margin")
```


traits.spe <- select(traits.sp, -sd)
traits.spe <- pivot_wider(traits.spe, names_from = Trait, values_from = mean)

#only keep corrected traits (divided by body length)
#only keep one head trait (keeping width)

traits.spe <- select(traits.spe, 1, 3,5,8, 9, 11, 13, 14)
traits.spe <- rename(traits.spe, species = X.1)
traits.spe <- as.data.frame(traits.spe) %>% ungroup()
traits.spe$species <- gsub(" ", "", traits.spe$species)
row.names(traits.spe) <- traits.spe$species
traits.spe <- select(traits.spe, -species)

spe.bin <- as.data.frame(spe.bin)
row.names(spe.bin) <- spe.bin$uniID
spe.bin <- select(spe.bin, -uniID)





#for traits species should be rows, traits cols
traits.spe <- select(traits.sp, -sd)
traits.spe <- pivot_wider(traits.spe, names_from = Trait, values_from = mean)

#only keep corrected traits (divided by body length)
#only keep one head trait (keeping width)

traits.spe <- select(traits.spe, 1, 3,5, 9, 11, 13, 14)
traits.spe <- rename(traits.spe, species = X.1)
row.names(traits.spe) <- traits.spe$species
traits.spe <- as.data.frame(traits.spe) %>% ungroup()
traits.spe <- select(traits.spe, -species)
#I want veg height but there are NAs ugh

env4 <- select(cov, uniID, 2, 3, 4, 5,12, 13,15,  17, 23, 29, 34)
env4$Microsite <- as.factor(env4$Microsite)

naveg <- env4[is.na(env4$veg..height),]

#make a vector of NA values
naveg <- naveg$uniID
`%!in%` <- Negate(`%in%`)

#remove those NA rows from the species table
spe <- spe[spe$uniID %!in% naveg, ]
env4 <- env4[env4$uniID %!in% naveg, ]

row.names(env4) <- env4$uniID
env4 <- env4[,-1]





m0  <- capscale(cwm.comm ~ 1+ Condition(Site), data = cwm.env, dist = "bray")



m1 <- capscale(cwm.comm ~ Prec + esi +mean.cover +dry.veg.percent+  + month+ Condition(Site), data = cwm.env, dist = "bray")


m1 <- capscale(cwm.comm ~ arid + esi  +mean.cover +dry.veg.percent+  Microsite  + Condition(Site), data = cwm.env, dist = "bray")

RsquareAdj(m1)$r.squared
m1
plot.cca(m1)


cor.test(cwm.env$esi, cwm.env$arid)
plot(cwm.env$esi, cwm.env$arid)

ggplot(cwm.env, aes(esi, arid)) + geom_smooth()

plot(m1, type = "n", scaling = "sites")
text(m1, dis = "cn", scaling = "sites")
text(m1, dis = "sp", scaling = "sites", col = "red")

autoplot(m1, scaling = 1)

summary(m1)
anova(m1)
anova.cca(m1, by = "terms", strata = cwm.env$Site)
anova.cca(m1, by = "margin", strata = cwm.env$Site)
anova.cca(m1, by = "axis", strata = cwm.env$Site)
ordiplot(m1)
arrows(m1, "species")
pl <- ordiplot(m1, type = "none")
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "species", col="blue", cex=0.9, arrows = TRUE)
arrows(p1, "sites")





env4$Site <- as.factor(env4$Site)
evn4C <- select(env4, 2, 3, 7, 8, 9, 11)
#evn4C <- select(env4, 2:11)
#evn4C <- select(evn4C, 1, 3)


#shrub and open are really close
#gotta add shrub type

library(mvabund)
fit <- traitglm(spe.bin, evn4C, traits.spe, family = "binomial", method = "manyglm")
fit$fourth #notice LASSO penalty has shrunk many interactions to zero

#summary(fit, blcok = env4)

anova(fit, block = env4$Site, resamp = "case", nBoot = 99)
#signficant trait*env but takes 36 minutes to run

library(lattice)
a        = max( abs(fit$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(fit$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

library(jtools)
summ(fit)
#test again for presence/absence
#ant abundance may be more related to colony size which is a trait




four.comb.ants.adj <- fourthcorner(env4, spe.bin,
                                   traits.spe, modeltype = 6, p.adjust.method.G = "none",
                                   p.adjust.method.D = "fdr", nrepet = 999)

four.comb.ants.adj 
plot(four.comb.ants.adj, alpha = 0.05, stat = "D2")
#why are the p adj so different for a bivariate correlation? 
#not sure


fit <- traitglm(spe.bin, env4, traits.spe,family = "binomial")
#summary(fit)
fit$fourth #notice LASSO penalty has shrunk many interactions to zero

a        = max( abs(fit$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(fit$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)


## Ant richness

```{r}
m1 <- glmmTMB(Species ~ arid +mean.cover +dry.veg.percent + Microsite  + shrub.site + month + (1|Site) , family = "poisson", data = cov)

m2 <- glmmTMB(Species ~ arid +dry.veg.percent + var.cover +Microsite  + shrub.site + month + (1|Site), family = "poisson", data = cov)

m3 <- glmmTMB(Species ~ arid +mean.cover +dry.veg.percent + var.cover +Microsite  + shrub.site + month + esi + (1|Site) , family = "poisson", data = cov)

m4 <- glmmTMB(Species ~ arid +mean.cover +dry.veg.percent  +Microsite  + shrub.site + month + esi + (1|Site) , family = "poisson", data = cov)



m5 <- glmmTMB(Species ~ mean.cover +dry.veg.percent  +Microsite  + shrub.site +  esi + arid + (1|Site) , family = "poisson", data = cov)
summary(m5)
AIC(m1,m2, m3, m4, m5)

ggplot(cov, aes(arid, Species)) + geom_smooth()



cor.test(cov$mean.cover, cov$var.cover)

summary(m2)
library(performance)

check_collinearity(m5)
m2 <- glmmTMB(Species ~ dry.veg.percent + Microsite  + mean.cover +  esi + mean.height +  shrub.site + arid + month + (1|Site) , family = "nbinom2", data = cov)
summary(m2)

AIC(m1, m2)
```





## Beta-diversity taxonomic






spe.bin <- as.data.frame(spe.bin)
row.names(spe.bin) <- spe.bin$uniID
spe.bin <- select(spe.bin, -uniID)

spe.bin.core <- betapart.core(spe.bin)
spe.bin.multi <- beta.multi(spe.bin)
spe.bin.samp <- beta.sample(spe.bin.core, sites = 9, samples = 100)


dist <- beta.pair(spe.bin, index.family = "sorensen")

bd <- vegan::betadisper(dist[[3]], cov$shrub.site)

boxplot(bd)
anova(bd)
# difference in dispersion between the three
TukeyHSD(bd)
# contrasts not significant. ok.



bd <- vegan::betadisper(dist[[3]], cov$Microsite)

boxplot(bd)
anova(bd)
TukeyHSD(bd)






library(FD)
library(vegan)
gow <- gowdis(env4)
euc <- dist(env4, "euclidean")
gow <- as.matrix(gow)
sor <- dist[[3]]
sor <- as.matrix(sor)

sor.mean <- rowMeans(sor, na.rm = TRUE)
sim <- as.matrix(dist[[1]])
sim.mean <- rowMeans(sim, na.rm = TRUE)
sne <- as.matrix(dist[[2]])
sne.mean <- rowMeans(sne, na.rm = TRUE)

envbeta <- cbind(env4, sor.mean, sim.mean, sne.mean)


mb <- glmmTMB(sim.mean ~  arid  + esi + mean.cover +dry.veg.percent+  Microsite  + shrub.site + (1|Site), data = envbeta)
shapiro.test(resid(msor))
plot(resid(mb))
summary(mb)
check_collinearity(mb)

ggplot(envbeta, aes(arid, sor.mean)) + geom_smooth() 


msor <- glmmTMB(sor.mean ~  arid  + esi + mean.cover +dry.veg.percent+  Microsite  + shrub.site + (1|Site), data = envbeta)
summary(msor)
msor2 <- glmmTMB(sor.mean ~  arid  + esi + mean.cover +dry.veg.percent+  Microsite  + shrub.site , data = envbeta)

summary(msor2)

AIC(msor, msor2)

plot(resid(mb) ~ envbeta$sim.mean)
mb2 <- glmmTMB(sne.mean ~  arid   +mean.cover +dry.veg.percent+  Microsite  + shrub.site +(1|Site), data = envbeta)

size <- filter(CWM.pres, Trait == "Webers")
m1 <- glmmTMB(cwm ~ esi + (1|Site.x), data = size)
summary(m1)

m2 <- glmmTMB(cwm ~ arid + (1|Site.x), data = size)
summary(m2)

ggplot(size, aes(arid, cwm)) + geom_smooth(method = lm)


summary(mb)
summary(mb2)


envlong <- pivot_longer(envbeta, 11:13, names_to = "beta", values_to = "mean")

envlong <- pivot_longer(envlong, c(1:6, 8, 9), names_to = "env", values_to = "value")

ggplot(envlong, aes(value,mean , color = beta)) +
  stat_smooth(method = "lm") + facet_wrap(~env, scales = "free")



euc.dry <-dist(env4$dry.veg.percent, "euclidean")
euc.Prec <- dist(env4$Prec, "euclidean")
euc.esi <- dist(env4$esi, "euclidean")
euc.vegheight <- dist(env4$veg..height, "euclidean")
euc.temp <- dist(env4$Temp, "euclidean")
euc.siteveg <- dist(env4$mean.bare, "euclidean")
euc.arid <- dist(env4$arid, "euclidean")

mantel(euc.dry, sor, method = "spearman", permutations = 999, na.rm = TRUE)


mantel(euc.Prec, sor, method = "spearman", permutations = 99, na.rm = TRUE)

mantel(euc.Prec, gow, method = "spearman", permutations = 99, na.rm = TRUE)

mantel(euc.vegheight, sor, method = "spearman", permutations = 99, na.rm = TRUE)

mantel(euc.temp, sor, method = "spearman", permutations = 99, na.rm = TRUE)

mantel(euc.siteveg, sor, method = "spearman", permutations = 99, na.rm = TRUE)

mantel(euc.arid, sor, method = "spearman", permutations = 99, na.rm = TRUE)
mantel(euc.esi, sor, method = "spearman", permutations = 99, na.rm = TRUE)

n <- mantel(euc.arid, sim, method = "spearman", permutations = 99, na.rm = TRUE)
mantel(euc.esi, sim, method = "spearman", permutations = 99, na.rm = TRUE)

plot(n)

# dry.veg. percent matters
# site level differences do

#yes, beta diversity (species composition) is correlated with the environment

## Beta-diversity functional

## db-RDA





# Sorenson

However, in this study we used
the C-score to evaluate co-occurrence between each pair
of species as C-score = (a − c) (b − c), where a is the total
number of occurrences of one species, b is the total number
of occurrences of a second species, and c is the number of
sites that contain both species




Contrarily, the Sørensen index
(Dice 1945) measures the mean number of shared sites (i.e.,
                                                      grid cells in our study) between a species pair, calculated
as S𝜙r = 2c∕(2c + a + b), where a, b and c are defined as
above. 




# Ant-habitat associations
## Indicator species analysis

- microsite
- shrub type at site
- study site

# Beta-diversity







## Functional richness



numberReps <- 100

#Lets create a matrix to store results from each iteration (one column per iteration)
resultsRandom4 <- matrix(NA, nrow = nrow(wide.pop), ncol = numberReps,
                         dimnames = list(rownames(wide.pop),
                                         paste0("Sim.", 1:numberReps)))
for(rep in 1:numberReps){
  randomizedFDis <- randomizeMatrix(samp = wide.pop, null.model = "independentswap")
  simFDis <- dbFD(trait.pop, randomizedFDis, w.abun = TRUE)$Q
  resultsRandom4[, rep] <- simFDis
}

obsFDis <- dbFD(trait.pop, wide.pop, w.abun = TRUE)$CWM
meanNull3 <- rowMeans(resultsRandom4)
ES3 <- obsFDis - meanNull3
sdNull3 <- apply(resultsRandom4, 1, sd)
SES3 <- ES3 / sdNull3
SES_rich <- data.frame(SES3)
SES_rich$sites <- row.names(SES_rich)
#join site-level environmental measurements to data

SES_rich <- rename(SES_rich, SES_ric = SES3)
SES <- left_join(SES_rich, SES, by = "sites")




mean(SES$SES_ric, na.rm = TRUE)

