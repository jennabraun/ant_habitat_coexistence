#cleaning

library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)


veg <- read.csv("raw data/vegetation.csv")

#loading the cleaned pitfall trap data (well, almost clean)

ants <- read.csv("D:/Research/Projects/prey_item_survey/Clean Data/arth_long.csv")
ants <- filter(ants, family == "Formicidae")
unique(ants$microsite)
unique(ants$site)
unique(ants$month)
unique(ants$rep)
#need to filter out NAs in these groups - I must have left them blank
#filter here so they don't count as zeroes

ants <- filter(ants, microsite != "NA" & microsite != "nolabel")
ants <- filter(ants, month != "NA")
ants <- filter(ants, site != "NA")
ants <- filter(ants, rep != "NA")


#I also want to filter out blank IDs & alates, but at this stage so I can sum easily
ants <- filter(ants, highest.rtu != "NA")
ants <- filter(ants, highest.rtu != "alate")

#need to extract a list of trap ID that counts insects but no ants
#making a unique identifier to facilitate joining tables
ants$uniID <- paste(ants$site, ants$month, ants$microsite, ants$rep)
ants$uniID <- gsub(" ","", ants$uniID)
reps <- unique(ants$uniID)


unique(ants$highest.rtu)
#less typos, not perfect yet

#ants <- filter(ants, site == "PaPl" | site == "CaS" | site == "Lok" | site == "SemiT" | site == "SiCr")
ants.ag <- dplyr::select(ants, uniID, highest.rtu, quantity)
#collapse and sum ants within species within rep
ants.ag <- ants.ag %>% group_by(uniID, highest.rtu) %>% summarise(quantity = sum(quantity)) 

#check to make sure everything adds correctin
sum(ants.ag$quantity, na.rm = TRUE)

#make a wide dataframe to calculate abundance and richness in vegan

wide <- ants.ag %>% spread(highest.rtu, quantity)
#replace NA with zero
wide[is.na(wide)] <- 0

#I made clean metadata with identifiers and the correct microsite in the prey item project already
#just bring in that csv and drop the arth data
metadata <- read.csv("D:/Research/Projects/prey_item_survey/Clean Data/arth_clean.csv")

#metadata <- metadata %>% filter(uniID %in% reps)

row.names(metadata) <- metadata$uniID
#weird idnames

#I need to extract the zeros from metadata and add to the bottom of wide

zeroes <- anti_join(metadata, wide, by = "uniID")
#90 without ants is realistic? I think I still need to remove the damaged pits, maybe not though

zeroes <- select(zeroes, "uniID")
wide <- bind_rows(wide, zeroes)
wide[is.na(wide)] <- 0
metadata <- filter(metadata, uniID != "NA")

metadata <- metadata[match(wide$uniID, metadata$uniID),]

all.equal(metadata$uniID, wide$uniID)

row.names(wide) <- wide$uniID
wide <- wide %>% ungroup(uniID) %>% select(-uniID)

metadata$abun <- apply(wide, 1, sum)
#check for total
sum(metadata$abun)
H <- diversity(wide)
simp <- diversity(wide, "simpson")
S <- specnumber(wide)
J <- H/log(S)
metadata$H <- H
metadata$Simpson <- simp
metadata$Species <- S
metadata$Even <- J


write.csv(metadata, "Clean Data/ants_clean.csv")
write.csv(wide, "Clean Data/wide_clean.csv")
write.csv(veg, "Clean Data/veg_clean.csv")
