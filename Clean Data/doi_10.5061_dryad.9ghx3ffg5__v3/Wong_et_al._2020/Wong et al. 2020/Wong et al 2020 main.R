library(sppairs)
library(Rmisc)
library(fBasics)
library(vegan)

source(file.choose()) #read species_summary2_revised.R
source(file.choose()) #read null_species_score_revised.R
source(file.choose()) #read OR_ratio2.R

data <- read.csv(file.choose()) #read All_abundance.csv
sites <- as.factor(data$sites)
invaded <- data$treatment
ant <- data[,c(-1,-2,-32),]

quantile(ant$Solenopsis_invicta[ant$Solenopsis_invicta > 0],c(.25,.50,.75)) #62,159,315 = 25,50,75th percentile
#ant$Solenopsis_invicta[ant$Solenopsis_invicta < 315] <- 0 #use this to convert number below threshold as zero

ant <- decostand(ant,"pa") #change it to presence and absence

### generate observed network here

OR.df <-spaa(ant,method="or.asymmetric2",asymmetric=T)
OR.df <- cbind(OR.df,log.OR=log(OR.df$value),sig=ifelse(test=(OR.df$value > 3 | OR.df$value < 1/3),yes="Sig",no="Insig")) #some studies have used observed OR > 3 / <1/3 as threshold for statistical significance. But we are not using it here (just FYI)

### generate null network here
set.seed(100)
library(vegan)

interaction.df <- list()
network_results<- species_results <- list()

for (j in 1 : 1000) { #generate 1000 null networks
  message("loop ", j)
  
  sm <- list()
  
  for (i in 1: length(levels(sites))) {
    site_df <- ant[which(sites == levels(sites)[[i]]),]
    sm[[i]] <- as.data.frame(simulate(nullmodel(site_df, "quasiswap")))
  }
  final_sm <- do.call(rbind,sm)
  
  colnames(final_sm) <- colnames(ant)
  
  null.df<- spaa(final_sm,method="or.asymmetric2",asymmetric=T)
  null.df <- cbind(null.df,log.OR=log(null.df$value),sig=ifelse(test=(null.df$value > 3 | null.df$value < 1/3),yes="Sig",no="Insig"),j=j)
  
  interaction.df[[j]] <- null.df
}

### species level z-score
null.df.link <- do.call(rbind,interaction.df)

species_null <- null_species_score(null.df.link,ant,all=T) #This one is null network
species_df <- null_species_score(cbind(OR.df,j=1),ant,all=T) #This one is observed network

z_score <- (species_df$average.diff - species_null$average.diff)/species_null$sd.diff #calculate z-score
names(z_score) <- colnames(ant) #name the ants

### link level z-score

pairwise_data_null <- null.df.link
pairwise_data_obs <- OR.df

null_link_df <- summarySE(data=pairwise_data_null, measurevar="log.OR", groupvars=c("sp1","sp2"))
SES_link <- (pairwise_data_obs$log.OR -null_link_df$log.OR)/null_link_df$sd
SES_link_df <- data.frame(sp1=null_link_df$sp1,sp2=null_link_df$sp2,SES=SES_link)
sinv_link <- subset(SES_link_df,sp1=="Solenopsis_invicta") #select links relevant to fire ants only

### phylogenetic tree
library(phytools)
library(phylolm)
library(rr2)
library(picante)

source(file.choose()) #read genus_tree2.R

tree <- read.tree(file.choose()) #read the ant_tree_FBD.tre from the Economo study
specific_df <- data.frame(z_score)

tree_set <- genus_tree(tree,specific_df,n=1)

pair.phy.dist <- final.result <- list()
for (k in 1:length(tree_set)) {
  pair.phy.dist[[k]] <- cophenetic.phylo(tree_set[[k]]) #calculate pairwise distance
  pair.phy.dist[[k]] <- pair.phy.dist[[k]][order(rownames(pair.phy.dist[[k]])),order(colnames(pair.phy.dist[[k]]))] #reorder tthe data
}

###environmental preference analyses
env <- read.csv(file.choose(),row.names=1) #read env_data.csv

env.sp.LMC <- env.sp.MP <- list()

for (i in 1:ncol(ant)) { #The environmental preference will vary according to the abundance threshold
  presence.num <- which(ant[,i] == 1)  
  env.sp.LMC[[i]] <- colMeans(env[presence.num[presence.num %in% which(sites == "LMC")],]) #env conditions in LMC for each species
  env.sp.MP[[i]] <- colMeans(env[presence.num[presence.num %in% which(sites == "MP")],]) #env conditions in MP for each species
}

df <- data.frame(do.call(rbind,env.sp.LMC),do.call(rbind,env.sp.MP))
rownames(df) <- colnames(ant)
env.sp <- data.frame(ndvi=1:29,temp=1:29,g.cover=1:29)
env.sp$ndvi <- apply(data.frame(df[,1],df[,4]),1,mean,na.rm=T)
env.sp$temp <- apply(data.frame(df[,2],df[,5]),1,mean,na.rm=T)
env.sp$g.cover <- apply(data.frame(df[,3],df[,6]),1,mean,na.rm=T)

rownames(env.sp) <- colnames(ant)
colnames(env.sp) <- colnames(env)

sinv.env.dist <- list()

for (n in 1:ncol(env.sp)){
  temp_env<-env.sp[,n]
  names(temp_env) <- rownames(env.sp)
  sinv.dist <- as.data.frame(as.matrix(vegdist(temp_env,method="manhattan",diag=T,upper=T))) #manhattan distance is equivalent to pairwise absolute difference 
  sinv.env.dist[[n]] <-  sinv.dist[-which(rownames(sinv.dist)=="Solenopsis_invicta"),"Solenopsis_invicta"] #only comparison with s. invcita is needed
}

sinv.env.dist.df <- as.data.frame(do.call(cbind,sinv.env.dist))
colnames(sinv.env.dist.df) <- colnames(env.sp)
rownames(sinv.env.dist.df) <- rownames(env.sp[-which(rownames(sinv.dist)=="Solenopsis_invicta"),])

test_df <- SES_link_df
test_df <- subset(test_df,sp1=="Solenopsis_invicta")

fitNDVI <- lm(test_df$SES~sinv.env.dist.df$ndvi)
fitTemp <- lm(test_df$SES~sinv.env.dist.df$temp)
fitgcover <- lm(test_df$SES~sinv.env.dist.df$g.cover)
summary(fitNDVI)
summary(fitTemp)
summary(fitgcover)

# this is for robust regression
library(robustbase)
fitNDVI <- lmrob(test_df$SES~sinv.env.dist.df$ndvi,control=lmrob.control(maxit.scale=1000,k.max=1000))
fitTemp <- lmrob(test_df$SES~sinv.env.dist.df$temp,control=lmrob.control(maxit.scale=1000,k.max=1000))
fitgcover <- lmrob(test_df$SES~sinv.env.dist.df$g.cover,control=lmrob.control(maxit.scale=1000,k.max=1000))
summary(fitNDVI)
summary(fitTemp)
summary(fitgcover)

### phylogenetic regression

final.result0 <- list()
for (k in 1:length(tree_set)) {
  tempTree <- tree_set[[k]]
  phy_df <- data.frame(SES_link_df,pair.phy.dist[[k]][pair.phy.dist[[k]] != 0]) #combine pairwise occurrence and phylogenetic distance data
  phy_df <- subset(phy_df,sp1=="Solenopsis_invicta") #only extract links relevant to s.invicta
  colnames(phy_df) <- c("sp1","sp2","SES","phy.dist")
  fit0 <- lm(SES~phy.dist,data=phy_df) 
  fit0 <- lmrob(SES~phy.dist,data=phy_df,control=lmrob.control(maxit.scale=1000,k.max=1000)) #can be run using robust regression
  final.result0[[k]] <- c(summary(fit0)$coefficients[,1],summary(fit0)$coefficients[,4],R2=summary(fit0)$r.squared,adj.R2=summary(fit0)$adj.r.squared)  #final.result8[[k]] <- c(summary(fit8)$coefficients[,1],summary(fit8)$coefficients[,4],R2=summary(fit8)$r.squared,adj.R2=summary(fit8)$adj.r.squared)
}

phylo_result_df0 <-as.data.frame(do.call(rbind,final.result0))

### trait regression
library(interactions)
trait<- read.csv(file.choose(),row.names=1) #read mean_trait.csv
pdf.trait <- read.csv(file.choose(),row.names=1) #read pdf_trait.csv

sinv <- trait[which(rownames(trait)=="Solenopsis_invicta"),]
diff <- trait[-which(rownames(trait)=="Solenopsis_invicta"),] - sinv[rep(1,each=nrow(trait[-which(rownames(trait)=="Solenopsis_invicta"),])),] #calculate HD as trait of the species minus trait of S. invicta
colnames(diff) <- c("abs.size","abs.leg","abs.head","abs.mand","abs.pron","abs.scap","abs.eye")

overall_df <- data.frame(test_df,diff,pdf.trait)
overall_df[,-1:-3] <- scale(overall_df[,-1:-3]) #standardize all trait predictor

## linear regression (Note the p value is not corrected; can conduct Bonferroni corrections by multiplying all by 6, which is the number of models being built)
fit1<-lm(SES~abs.pron*pdf.pron,data=overall_df)
summary(fit1)
fit2<-lm(SES~abs.head*pdf.head,data=overall_df)
summary(fit2)
fit3<-lm(SES~abs.size*pdf.size,data=overall_df)
summary(fit3)
fit4<-lm(SES~abs.leg*pdf.leg,data=overall_df)
summary(fit4)
fit5<-lm(SES~abs.eye*pdf.eye,data=overall_df)
summary(fit5)
fit6<-lm(SES~abs.scap*pdf.scap, data=overall_df)
summary(fit6)

# can also be run using robust regression
#fit1<-lmrob(SES~abs.pron*pdf.pron,data=overall_df,control=lmrob.control(maxit.scale=1000,k.max=1000))
#fit2<-lmrob(SES~abs.head*pdf.head,data=overall_df,control=lmrob.control(maxit.scale=1000))
#fit3<-lmrob(SES~abs.size*pdf.size,data=overall_df,control=lmrob.control(maxit.scale=1000))
#fit4<-lmrob(SES~abs.leg*pdf.leg,data=overall_df,control=lmrob.control(maxit.scale=1000))
#fit5<-lmrob(SES~abs.eye*pdf.eye,data=overall_df,control=lmrob.control(maxit.scale=1000))
#fit6 <- lmrob(SES~abs.scap*pdf.scap, data=overall_df,control=lmrob.control(maxit.scale=1000))

overall_df <- data.frame(test_df,diff,pdf.trait) #rerun the model without standardization
fit1<-lm(SES~abs.pron*pdf.pron,data=overall_df) 
bounds <- johnson_neyman(fit1,abs.pron,pdf.pron,plot=F,control.fdr = T)
bounds

##test non-linearity 
fit1a<-lm(SES~poly(abs.pron,2)*poly(pdf.pron,2),data=overall_df)
summary(fit1a)
fit1b<-lm(SES~poly(abs.pron,2)*poly(pdf.pron,1),data=overall_df)
summary(fit1b)
fit1c<-lm(SES~poly(abs.pron,1)*poly(pdf.pron,2),data=overall_df)
summary(fit1c)
fit1d<-lm(SES~poly(abs.pron,1)*poly(pdf.pron,1),data=overall_df)
summary(fit1d)
fit1e<-lm(SES~poly(abs.pron,1)+poly(pdf.pron,1),data=overall_df)
summary(fit1e)
fit1f<-lm(SES~poly(abs.pron,2)+poly(pdf.pron,1),data=overall_df)
summary(fit1f)
fit1g<-lm(SES~poly(abs.pron,1)+poly(pdf.pron,2),data=overall_df)
summary(fit1g)

library(MuMIn)
AICc(fit1a,fit1b,fit1c,fit1d,fit1e,fit1f,fit1g) #1d is the best