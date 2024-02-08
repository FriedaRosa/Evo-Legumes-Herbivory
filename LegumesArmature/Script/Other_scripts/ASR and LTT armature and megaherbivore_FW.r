rm(list=ls())
library(geiger)
library(phytools)

#file.choose()
setwd("~/iDiv Dropbox/Renske Onstein/My Mac (idivmac32.local)/Documents/Project Legume traits/R_database")
traitData <- read.csv("Matrix_all_traits.csv", stringsAsFactors = FALSE, header = TRUE)
str(traitData)
summary(traitData)
dd <-traitData  %>% filter(!is.na(Armatures)) #select only those without NA in Armature (this is also important to not have NA in the tree)

# Load phylogeny
mimosoidPhylo <- read.tree("/Users/ro68gizi/iDiv Dropbox/Renske Onstein/My Mac (idivmac32.local)/Documents/Project Legume traits/Erik tree/Mimo_metachronogram_mmc.tre")

#adjust tip names
new_tip_names<-read.csv("/Users/ro68gizi/iDiv Dropbox/Renske Onstein/My Mac (idivmac32.local)/Documents/Project Legume traits/Erik tree/replace_tip_labels.csv")
tree <- mimosoidPhylo

# Match species names in data with those in tree:
tree_tips <- unique(tree$tip.label) #name vector with tip-labels in tree
setdiff(tree_tips, new_tip_names$Tiplabels_metatree) #identical

dd2 <- merge(dd, new_tip_names, by.x="Taxon", by.y="Taxanames_dataset") # merge tree-tip-names with trait data
dd2$Taxon <- NULL # remove old names


# Identify species that differ between the phylogeny and trait dataset
setdiff(tree$tip.label, dd2$Tiplabels_metatree)
setdiff(dd2$Tiplabels_metatree, tree$tip.label)

# remove duplicated species and only keep unique entries
dd3 <- dd2 %>% 
  group_by(Tiplabels_metatree) %>%
  mutate(duplicate = n()) %>% # count number in each group
  filter(duplicate == 1) %>% # select only unique records
  select(-duplicate) # remove group count column

dd3 <- as.data.frame(dd3)

# Exclude those taxa from the phylogeny and the trait dataset
intersectTaxa <- intersect(dd3$Tiplabels_metatree, tree$tip.label)
traitDataSubset <- subset(dd3, Tiplabels_metatree %in% intersectTaxa)
mimoPhyloSubset <- drop.tip(tree,
                            tip = tree$tip.label[!tree$tip.label %in% intersectTaxa])



# match data and tree ===================================
row.names(dd3) <- dd3$Tiplabels_metatree
setdiff(tree$tip.label, dd3$Tiplabels_metatree)

# Drop these species from the tree, because we need matching data between tree and traits
tree2 <- drop.tip(tree, setdiff(tree$tip.label, dd3$Tiplabels_metatree))
tree2 <- ladderize(tree2)

matches <- match(dd3$Tiplabels_metatree, tree2$tip.label, nomatch = 0)
data2 <- subset(dd3, matches != 0)
data2 <- ReorderData(tree2, data2, taxa.names = "row names")
data2$Tiplabels_metatree <- NULL
dd3 <- data2
str(dd3) # this data frame has the same species names as those in the tree and is ordered by the same order as the tip labels in the tree

# Select Armature trait (and make a vector from it)
armature <- numeric(length = nrow(dd3))
names(armature) <- rownames(dd3)
unique(dd3$Armatures)
armature[dd3[,"Armatures"]=="0"] <- "absent"
armature[dd3[,"Armatures"]==1] <- "present"
armature <- as.factor(armature)
levels(armature)

armature <- ReorderData(tree2, armature, taxa.names="names") # ordered trait vector for models (not dataframe)
str(armature) # make sure it's a factor

tree <- tree2
trait <- armature
###########################################################################################
# Stochastic character mapping

#Mk1 all rates equal (do you remember this model? What does it do?)
equal <- fitDiscrete(tree, trait, model = "ER")
#Mk all rates different (do you remember this model? What does it do?)
ard <- fitDiscrete(tree, trait, model = "ARD")

# Calculate AIC weights again
aic.discrete <- setNames(c(equal$opt$aic, ard$opt$aic), c("equal", "different"))
weights <- aicw(aic.discrete)
weights

# The ER model fits best, and we use this to infer ancestral states
ans<-ace(trait, tree, type = "discrete", marginal = FALSE, model = "ER")
plot(tree, label.offset = 1)
co <- c("blue", "yellow") # armature present (blue), absent (yellow)
tiplabels(pch = 22, bg = co[as.numeric(trait)], cex = 2, adj = 1)
nodelabels(pie = ans$lik.anc, piecol = co, cex = 0.3)

# simmap
mtree <- make.simmap(tree, trait, model = "ER") #under ARD or ER model
# in case you know which state the root has, you can add this information with the additional argument: e.g. "pi=c(1,0))" (meaning the root has 100% probability for first trait state (alphabetically) and 0 % probability for the second trait state)

cols <- setNames(c("black", "gold"), sort(unique(trait)))
plotSimmap(mtree, cols, pts = FALSE, lwd = 2)
add.simmap.legend(colors = cols, vertical = FALSE, prompt = FALSE, x = 0, y = 24)

#A single stochastic character map does not mean a whole lot in isolation; 
#however in aggregate, we can do quite a bit. For instance, we can estimate 
#the number of changes of each type, the proportion of time spent in each 
#state, and the posterior probabilities that each internal node is in each 
#state, under our model. For example:

#For 100 simulation (nsim=100)
mtrees <- make.simmap(tree, trait, model = "ER", nsim = 100)
XX2 <- describe.simmap(mtrees, plot = FALSE)

#Can take a while so save your workspace afterwards:
save.image("XXX")

## now let's plot a random map, and overlay the posterior probabilities
plotSimmap(mtrees[[1]], cols, lwd = 1, pts = F, setEnv = T)
#plotSimmap(mtree, cols, pts = FALSE, lwd = 2)

#Change colours
cols2 <- setNames(c("white", "gold"), sort(unique(trait)))
#nodelabels(pie = XX2$ace, piecol = cols2, cex = 0.3) #this can be done if you also did a simple ARD model in e.g. ace
add.simmap.legend(colors = cols, vertical = FALSE, prompt = FALSE, x = 0, y = 24)

#We can also use stochastic mapping to plot the posterior probability that the edges & nodes of the 
#tree are in a binary state. This function is somewhat primitive, so we first 
#have to convert our mapped edge states to 0 & 1.

## convert mapped states
stateb <- mergeMappedStates(mtrees, c("0"), "0")
stateb <- mergeMappedStates(stateb, "1", "1")
## now plot density map
XX <- densityMap(stateb, lwd = 2)
obj<-XX

## what is the length of the current color ramp?
n<-length(obj$cols)

## change to black -> gold
obj$cols[1:n]<-colorRampPalette(c("black","gold"), space="Lab")(n)
plot(obj, lwd=1.5)

## change to grey -> gold
obj$cols[1:n]<-colorRampPalette(c("grey","gold"), space="Lab")(n)
plot(obj, lwd=2, show.tip.label=FALSE)
axisPhylo()

save.image("XXX")



###########################################################################################

#ASR and diversity estimates through time

library(ape)
library(geiger)
library(phytools)

setwd("~/iDiv Dropbox/Renske Onstein/My Mac (idivmac32.local)/Documents/iDiv/students/PhD/Ferreira_Rachel/analysis")

#You can use this function to recreate many of the figures from our 2010 Evolution paper. 
#The following script gives an example for the plot that illustrates how lineage diversity estimates for nodes vary over time, 
#as anoles colonize different islands and radiate (different islands have different colors). This script assumes we've already read in a time-calibrated tree (GA_Anolis), 
#and that we have a vector indicating the geographic location of each species (geography).

## First estimate ancestral area for each node, using unequal rate model (ARD) or equal rate moder (ER)
recon <- rerootingMethod(tree= tree,x=trait,model="ER")
for.piecharts <- recon$marginal.anc
write.csv(for.piecharts, file="for.piecharts_ER_mimosoid.csv")

#plot the marginal ancestral states
#0=no armature, 1=armature
col <- c("grey","gold")

#ARD palms
tree3<-ladderize (tree_mimosoid)
plot(tree3, cex=0.3)
nodelabels(pie=for.piecharts, piecol=col, cex=.3)
legend("bottomleft", c("no armature", "armature"), fill=col)

#Remove negative values in the ARD model
for.piecharts3<-for.piecharts
for.piecharts3[for.piecharts3<0]<-0
for.piecharts4<-for.piecharts2
for.piecharts4[for.piecharts4<0]<-0

## Calculate ages of all nodes in the trees.
time.vector2 <- max(branching.times(tree_mimosoid))-branching.times(tree_mimosoid)
time.vector2<-branching.times(tree_mimosoid)

## Estimate lineage diversities of all nodes in the tree under model ARD and model ER (default)
diversity.vector <- estDiversity(tree_mimosoid,x=armature,method="asr", model="ER")
write.csv(diversity.vector, file="diversity.vector_ER_primates.csv")

save.image("mimosoid est diversity.RData")

## Make a vector of colors for the areas (in the example, we've got 4 islands).
piecolors <- c("#004165", "#eaab00")
piecolors <-col <- c("grey","gold")

## Finally, make a plot of lineage diversity estimates over time for all nodes.
require(mapplots)
plotNodePie = function(xvals,yvals,pie,pie.colors,radius=radius,...){
  plotPie = function(temp,pie.colors){
    add.pie(x=temp["xvals"],y=temp["yvals"],z=temp[-c(1:2)],radius=radius,col=pie.colors,labels=NA,init.angle=0)
  }
  plot(xvals,yvals,type="n",...)
  temp.matrix = cbind(xvals,yvals,pie)
  apply(temp.matrix,MARGIN=1,plotPie,pie.colors=pie.colors)
  return("done")
}

par(mfrow=c(2,1))

pdf("plot.pdf", width = 16, height = 7)
dev.off()

#Mimosoid armature
plotNodePie(xvals= time.vector2,yvals= diversity.vector2,pie=for.piecharts2,pie.colors=c("grey", "red"),xlab="Time From Root",ylab="Lineage Diversity Estimate",radius=25, xlim=c(120, 0))
legend("topleft",c("no armature", "armature"),fill=c("grey", "red"))

