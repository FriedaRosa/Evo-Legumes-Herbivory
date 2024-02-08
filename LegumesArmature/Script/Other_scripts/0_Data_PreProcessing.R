rm(list=ls())
gc()
setwd("../")

## Output path
out.path <- "output/"


## Libraries 
library(geiger)
library(phytools)
library(dplyr)
library(evobiR)
library(phangorn)


## Script to create final MCC tree and data for A) Mammals and B) Legumes =====  #

# [A] Mammals =================================================================  #

## 1.] Compute MCC tree from 1000 trees

p.mammals.1000trees <- "input/Mammals/Raw/mammal_1000.nex" # 1000 trees from posterior distribution of mammal phylogeny
tree_mammals <- read.nexus(p.mammals.1000trees) # Raw Data

MCC_mammals <- mcc(tree_mammals, tree = TRUE, part = NULL, rooted = TRUE) # Data Processing

# write.tree(MCC_mammals, paste(p.mammals.MCC)) # Procesed Data (MCC tree)

rm(tree_mammals) # remove raw data








# 2.] Create Model dataframe with species and traits
p.mammals.traits <- "input/Mammals/Raw/Trait_data.csv" # all traits for all mammals
mammal_traitData <- read.csv(p.mammals.traits, sep=",")

##make dataframe with species, mass  trait
mammal_traitData2 <- mammal_traitData[,c(1,6,12,20)]
mammal_subset_df <- mammal_traitData2[mammal_traitData2[, 2] == 1, ]

# Initialize a new column filled with zeros
mammal_subset_df$Condition_Met <- 0 

# Update this column to 1 for animals that meet both conditions (mass >= 10000 and diet >= 95)
mammal_subset_df$Condition_Met[mammal_subset_df$Mass.g >= 10000 & mammal_subset_df$Diet.Plant >= 95] <- 1

mammal_traitData3 <- mammal_subset_df[,c(1,5)]
str(mammal_traitData3)

# write.csv(mammal_traitData, file = paste(p.mammals.traits.short), row.names = FALSE)
rm(mammal_subset_df, mammal_traitData,mammal_traitData2)

mammal_traitData <- mammal_traitData3
rm(mammal_traitData3)







# Exclude those taxa from the phylogeny and the trait dataset
mammal_intersectTaxa <- intersect(mammal_traitData$Binomial.1.2, MCC_mammals$tip.label)

# Drop species from the trait data
mammal_traitDataSubset <- subset(mammal_traitData, Binomial.1.2 %in% mammal_intersectTaxa)
mammal_traitDataSubset2 <- data.frame(row.names = mammal_traitDataSubset$Binomial.1.2, 
                                      Condition_Met = mammal_traitDataSubset$Condition_Met)
# write.csv(mammal_traitDataSubset2, paste(p.mammals.traits.MCC))

# Drop tips from the tree
mammal_PhyloSubset <- drop.tip(MCC_mammals,
                               tip = MCC_mammals$tip.label[!MCC_mammals$tip.label %in% mammal_intersectTaxa])
# write.tree(mammal_PhyloSubset, paste(p.mammals.MCC_drop))
rm(mammal_traitDataSubset, mammal_intersectTaxa)








# Select Mass trait (and make a vector from it)
Mass <- numeric(length = nrow(mammal_traitDataSubset2))

# Create the 'Mass' factor based on the 'Condition_Met' column
Mass <- ifelse(mammal_traitDataSubset2$Condition_Met == 0, "other", "megaherbivore")

# Convert to a factor
Mass <- as.factor(Mass)
table(Mass)

# Reorder data if needed (assuming ReorderData and tree2 are defined)
Mass <- ReorderData(mammal_PhyloSubset, Mass, taxa.names="names")
names(Mass) <- rownames(mammal_traitDataSubset2)

# ordered trait vector for models (not dataframe)
str(Mass) # make sure it's a factor

saveRDS(Mass, "input/Mammals/Processed/Mass.rds")










# [B] Legumes =================================================================  #

p.legumes.traits <- "input/Legumes/Raw/Matrix_all_traits.csv" # all traits for all legumes
p.legumes.traits.corrected <- "input/Legumes/Processed/Matrix_all_traits_tipsmatched.csv"







## Match species names in data with tree
legumes_traitData <- read.csv(p.legumes.traits, stringsAsFactors = FALSE, header = TRUE)

legumes_traitData <- legumes_traitData %>% 
  # select(Taxon, Armatures) %>% 
  distinct(.)
table(legumes_traitData$Armatures)

legumes_dd <- legumes_traitData  %>% 
  filter(!is.na(Armatures)) #select only those without NA in Armature (this is also important to not have NA in the tree)

# Load phylogeny
p.legumes.MCC <- "input/Legumes/Raw/Mimo_metachronogram_mmc.tre" # Legumes MCC tree
MCC_legumes <- read.tree(p.legumes.MCC)

## adjust tip names
#  Tip labels to match data and phylo
p.legumes.new.tips <- "input/Legumes/Raw/replace_tip_labels.csv"
new_tip_names <- read.csv(p.legumes.new.tips)

# Match species names in data with those in tree:
tree_tips <- unique(MCC_legumes$tip.label) #name vector with tip-labels in tree

# Double-check if names match:
setdiff(tree_tips, new_tip_names$Tiplabels_metatree) #identical
setdiff(new_tip_names$Tiplabels_metatree, tree_tips) #identical


# merge tree-tip-names with trait data
legumes_dd2 <- merge(legumes_dd, new_tip_names, by.x="Taxon", by.y="Taxanames_dataset") 
# write.csv(legumes_dd2, "input/Matrix_all_traits_tipsmatched.csv")

legumes_dd2$Taxon <- NULL # remove old names

legume_traits <- legumes_dd2 %>% select(Tiplabels_metatree, Armatures)
# write.csv(legume_traits, paste0(out.path, "CSV", "legume_traits_short.csv"))

rm(legumes_traitData, legumes_dd, new_tip_names, tree_tips, legumes_dd2)








# Exclude those taxa from the phylogeny and the trait dataset
mimo_intersectTaxa <- intersect(legume_traits$Tiplabels_metatree, MCC_legumes$tip.label)

mimo_traitDataSubset <- subset(legume_traits, Tiplabels_metatree %in% mimo_intersectTaxa) %>% 
  select(Tiplabels_metatree, Armatures) %>% 
  distinct(.)
mimo_traitDataSubset2 <- data.frame(row.names = mimo_traitDataSubset$Binomial.1.2, 
                                    Armatures = mimo_traitDataSubset$Armatures)

mimo_PhyloSubset <- drop.tip(MCC_legumes,
                             tip = MCC_legumes$tip.label[!MCC_legumes$tip.label %in% mimo_intersectTaxa])

rm(mimo_traitDataSubset, mimo_intersectTaxa)








# Select Armatures trait (and make a vector from it)
Armature <- numeric(length = nrow(mimo_traitDataSubset2))

# Create the 'Armature' factor based on the 'Armatures' column
Armature <- ifelse(mimo_traitDataSubset2$Armatures == 0, "no_armature", "armature")

# Convert to a factor
Armature <- as.factor(Armature)
table(Armature)

# Reorder data if needed (assuming ReorderData and tree2 are defined)
Armature <- ReorderData(mimo_PhyloSubset, Armature, taxa.names="names")
names(Armature) <- rownames(mimo_traitDataSubset2)

# ordered trait vector for models (not dataframe)
str(Armature) # make sure it's a factor

saveRDS(Armature, "input/Legumes/Processed/Armature.rds")







### Summary ===============

# 1.] MCC Trees
MCC_legumes

MCC_mammals

## Pruned:
mimo_PhyloSubset 

mammal_PhyloSubset

# 2.] filtered trait data
legume_traits
mammal_traitData

# Pruned:
mimo_traitDataSubset2
table(mimo_traitDataSubset2)

mammal_traitDataSubset2
table(mammal_traitDataSubset2)

## Save to input/Processed folders
saveRDS(mimo_PhyloSubset, paste0("input/Legumes/Processed/Legumes_MCC_drop.rds"))
saveRDS(mimo_traitDataSubset2, paste0("input/Legumes/Processed/Legumes_traits_drop.rds"))

saveRDS(mammal_PhyloSubset, paste0("input/Mammals/Processed/Mammals_MCC_drop.rds"))
saveRDS(mammal_traitDataSubset2, paste0("input/Mammals/Processed/Mammals_traits_drop.rds"))


