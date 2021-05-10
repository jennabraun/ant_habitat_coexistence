#cleaning

library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)

veg <- read.csv("raw data/vegetation.csv")
str(veg)
#this is a check for typos in factor names
unique(veg$month)
unique(veg$Method)
unique(veg$Microsite)
unique(veg$Site)

#need to fill in the field microsite field - field micro is the labels on the vials and corresponds to
#the insect data, but Microsite is the actual microsite - techs labelled shrub/open at shrubless sites to aid organization
veg$field.micro
veg <- veg %>% mutate(field.micro = ifelse(nchar(veg$field.micro) ==0, paste(veg$Microsite, veg$Rep), field.micro))
veg$field.micro <- gsub(" ","", veg$field.micro)


#loading pitfall data from another folder, this is ongoing so pulling direct from excel
library(XLConnect)
ants <- readWorksheetFromFile("D:/Research/Projects/prey_item_survey/Raw Data/Arthropod_data.xlsx", "Pitfalls")

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

#analyze ants only
ants <- ants %>% filter(family == "Formicidae")
sum(ants$quantity, na.rm = TRUE)
unique(ants$genus)
ants$genus <- gsub(" ", "", ants$genus)
unique(ants$highest.rtu)


#ants <- filter(ants, site == "PaPl" | site == "CaS" | site == "Lok" | site == "SemiT" | site == "SiCr")
ants.ag <- dplyr::select(ants, uniID, genus, quantity)
#collapse and sum ants within species within rep
ants.ag <- ants.ag %>% group_by(uniID, genus) %>% summarise(quantity = sum(quantity)) 

#check to make sure everything adds correctin
sum(ants.ag$quantity, na.rm = TRUE)

#make a wide dataframe to calculate abundance and richness in vegan

wide <- ants.ag %>% spread(genus, quantity)
#replace NA with zero
wide[is.na(wide)] <- 0

#need to create matching unique identifier for the vegetation data and join

metadata <- veg
metadata$uniID <- paste(metadata$Site, metadata$month, metadata$field.micro)
metadata$uniID <- gsub(" ","", metadata$uniID)

#filter out unsampled pitfalls but keep the ones with zero ants

metadata <- metadata %>% filter(uniID %in% reps)

row.names(metadata) <- metadata$uniID
#weird idnames

#I need to extract the zeros from metadata and add to the bottom of wide

zeroes <- anti_join(metadata, wide, by = "uniID")
#150 without ants seems like a lot but ok
#ahh that was due to NAs 44 is much more realistic
str(zeroes)

zeroes <- select(zeroes, "uniID")
wide <- bind_rows(wide, zeroes)
wide[is.na(wide)] <- 0
metadata <- filter(metadata, uniID != "NA")
wide <- filter(wide, uniID != "NA")
#ok there are 4 rows that don't match for some unknown reasons

#weird <- anti_join(wide, metadata, by = "uniID")


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
