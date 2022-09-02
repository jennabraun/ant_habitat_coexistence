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
