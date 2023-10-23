rm(list=ls())
library(geiger)
library(phytools)
library(dplyr)
library(evobiR)
library(RColorBrewer)
library(tictoc)



# Read in the data and tree ============================

traitData <- read.csv("Matrix_all_traits.csv", stringsAsFactors = FALSE, header = TRUE)
tree <- read.tree("Mimo_metachronogram_mmc.tre")
tree <- as.phylo(tree)
new_tip_names<-read.csv("replace_tip_labels.csv")

# Subset data by trait 
armature <- traitData %>% select(Taxon, Armatures) %>% filter(!is.na(Taxon) & !is.na(Armatures)) %>% distinct(.) # 2 species were duplicated (Acacia baileyana and Chidlowiana sanguinea)

# Correct tip names in the tree with new labels
tree2 <- tree
tree2$tip.label <- new_tip_names[[2]][match(tree$tip.label, new_tip_names[[1]])]

# match data and tree ===================================
row.names(armature) <- armature$Taxon
sp_drop <- setdiff(tree2$tip.label, armature$Taxon)

# Drop these species from the tree, because we need matching data between tree and traits
tree3 <- drop.tip(tree2, sp_drop)
tree4 <- ladderize(tree3)

matches <- match(armature$Taxon, tree4$tip.label, nomatch = 0)
armature2 <- subset(armature, matches != 0)

setdiff(tree4$tip.label, armature2$Taxon) # 0
setdiff(armature2$Taxon, tree4$tip.label) # 0

tree <- tree4

# Make data for Mk models for Armature

armature2[armature2=="0"] <- "no_armature"
armature2[armature2=="1"] <- "armature"

trait <- unique(armature2)
trait <- ReorderData(tree, trait, taxa.names = "row names")

trait$Armatures <- as.factor(trait$Armatures)
trait$Taxon <- NULL

# Clean environment to free memory:
rm(tree2, tree3, tree4, traitData, new_tip_names, armature, matches, sp_drop)

# Save workspace in between =======================
save.image("Fig4Pollination/workspace_before_Mk_models_MCC.RData")
load("../Fig4Pollination/Fig4Pollination/workspace_after_Mk_models_MCC.RData")

# Mk models ============================
# estimate the rates of transitioning between ancestral states given a hypothesis (either equal rates(ER), or different rates(ARD))

# computing time = 4.67s
system.time(
  equal <- fitDiscrete(tree, trait, model = "ER"), gcFirst = T
  )

# computing time = 17.55s
system.time(
  ard <- fitDiscrete(tree, trait, model = "ARD")
  )

# Calculate AIC weights
aic.discrete <- setNames(c(equal$opt$aic, ard$opt$aic), c("equal", "different"))
weights <- aicw(aic.discrete)
weights # ARD model performs better for armature/no armature

# we need a slightly different format for Simmap:
trait2 <- as.factor(trait$Armatures)
names(trait2) <- rownames(trait)

# For 100 simulation (nsim=100) --> computing time = 167.91s
system.time(
  mtrees <- make.simmap(tree, trait2, model = "ARD", nsim = 100)
)

# 2.5s
system.time(
  pd_mtrees <- describe.simmap(mtrees, plot = FALSE)
)

# Save workspace in between =======================
save.image("Fig4Pollination/workspace_after_Mk_models_MCC.RData")
load(file = "Fig4Pollination/workspace_after_Mk_models_MCC.RData")

# We'd be interested in tracking the number of changes per unit of edge length rather than the total number of changes
# - because in all reconstructed phylogeny of extant taxa there is more edge length towards the tips of the tree than towards the root.
# we do this using the following argument: type = "rate"

# computing time = takes > 12h
tic()
  obj <- ctt(mtrees[1], type="rate")
toc()

###### Script from Pollination study: https://github.com/rubysaltbush/pollination-macroevolution/


# this function comes from https://github.com/rubysaltbush/pollination-macroevolution/blob/main/scripts/functions/transition_times.R

transition_times <- function(simmap){
  # below adapted from Liam Revells' phytools blog 
  # http://blog.phytools.org/2015/08/getting-timing-of-trait-changes-from.html
  # extracts raw transition times from a simmap (collapses multiple transitions
  # down into single transition events)
  # get tips and their states
  
  
  
  x <- phytools::getStates(simmap,"tips")
  levs<-sort(unique(c(getStates(tree,"tips"),
                      getStates(tree,"nodes"))))
  ct<-map.to.singleton(tree)
  
  
  # 
  H <- nodeHeight(ct)
  h<-c(0,max(H)-branching.times(ct),min(sapply(1:Ntip(ct),
                                               nodeheight,tree=ct)))
  ss<-setNames(as.factor(names(ct$edge.length)),
               ct$edge[,2])
  lineages<-matrix(0,length(h),length(levs),
                   dimnames=list(names(h),levs))
  lineages[1,getStates(tree,"nodes")[1]]<-1
  for(i in 2:length(h)){
    ii<-intersect(which(h[i]>H[,1]),which(h[i]<=H[,2]))
    lineages[i,]<-summary(ss[ii])
  }
  
  rm(x, states)

  ii<-order(h)
  times<-h[ii]
  

  lineages<-lineages[ii,]
  lineages<-cbind(lineages,total=rowSums(lineages))
  obj<-list(times=times,ltt=lineages)	
  if(gamma==FALSE){
    obj<-list(ltt=lineages,times=times,tree=tree)
    class(obj)<-"ltt.simmap"
  } else {
    gam<-gammatest(ltt(as.phylo(tree),plot=FALSE))
    obj<-list(ltt=lineages,times=times,gamma=gam$gamma,
              p=gam$p,tree=tree)
    class(obj)<-"ltt.simmap"
  }
  
}
  
  
  transition_times(mtrees)
  
  
  
  # then a list to store results in
  changes <- vector(mode="list", length = m*(m - 1))
  rm(m)
  # named by types of transitions
  names(changes) <- as.vector(ct[ii])
  rm(ct, ii)
  # then singling out maps where transitions happen (where there is more than 1 state)
  nc <- sapply(simmap$maps, length) - 1
  ind <- which(nc > 0)
  nc <- nc[ind]
  
  # getting the node heights (measure of time/branch lengths) across the tree
  H <- phytools::nodeHeights(simmap)
  maps <- simmap$maps[ind]
  # then looping through and calculating the node heights of each transition
  for(i in 1:length(maps)){
    for(j in 1:nc[i]){
      sc <- paste(names(maps[[i]])[j:(j + 1)], collapse = "->")
      h <- H[ind[i], 1] + cumsum(maps[[i]])[j]
      changes[[sc]] <- c(changes[[sc]], as.numeric(h))
    }
  }
  rm(nc, ind, h, H, i, j, sc, maps)
  # removing any nulls from list of changes and sorting small to large
  changes <- changes[!sapply(changes, is.null)]
  changes <- lapply(changes, sort, decreasing = FALSE)
  
  # now convert this changes list into nice data frame output
  output <- data.frame()
  for(i in 1:length(changes)){
    df <- dplyr::bind_cols(changes[i])
    df <- df %>%
      mutate(transition = colnames(df)) %>%
      rename(nodeheight = 1)
    output <- rbind(output, df)
  }
  
  # node heights are the height above the root, so time but inverse along the tree
  # to get time from node heights need to subtract from max height of tree
  output$time <- max(nodeHeights(simmap)) - output$nodeheight
  
  # get rid of nodeheight column
  output <- output[-1]
  
  # and return the output! to graph etc.
  output
}

## below adapted from other script:  https://github.com/rubysaltbush/pollination-macroevolution/blob/main/scripts/analysis/simmap.R

# apply function across list of multiple simulations
armature_transitions <- data.frame()
for(i in 1:length(mtrees)){
  temp <- cbind(i, transition_times(mtrees[[i]]))
  armature_transitions <- rbind(armature_transitions, temp)
}
rm(temp, i)

table(armature_transitions$transition)

# build new data frame with cumulative number of transitions
armature_trans_cumul <- data.frame()
for(n in 1:100){
  trans <- armature_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  armature_trans_cumul <- rbind(armature_trans_cumul, trans)
}
rm(n, trans)


##### copied from somewhere else in the code from the same script....:

# first need to know average number of transitions, rounded
trans_avg_length <- armature_trans_cumul %>%
  dplyr::group_by(transition, i) %>%
  dplyr::mutate(no_trans = max(trans_no)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(trans_no == no_trans) %>%
  dplyr::group_by(transition) %>%
  dplyr::mutate(avg_length = round(mean(no_trans), digits = 0)) %>%
  dplyr::select(transition, avg_length) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

# now knowing this, can rearrange data and average times across rows
avg_trans_times <- armature_trans_cumul %>%
  dplyr::group_by(transition, trans_no) %>%
  dplyr::mutate(avg_time = mean(time)) %>%
  dplyr::mutate(SE_time = sqrt(var(time) / length(time))) %>%
  dplyr::select(transition, trans_no, avg_time, SE_time) %>%
  dplyr::distinct()

# reduce avg_trans_times to trans_avg_length
avg_trans_times_a2na <- avg_trans_times %>%
  dplyr::filter(transition == "armature->no_armature") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
avg_trans_times_na2a <- avg_trans_times %>%
  dplyr::filter(transition == "no_armature->armature") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])
avg_trans_times <- rbind(avg_trans_times_a2na, avg_trans_times_na2a)
rm(avg_trans_times_a2na, avg_trans_times_na2a)

# export these results to csv in case I need them
readr::write_csv(avg_trans_times, "../Fig4Pollination/mean_transition_times_armature_MCC.csv")



#### Plotting ======================================

# install.packages("prettyGraphs") #to set transparancy of colors with "alpha"
library(prettyGraphs)

myColours = c("steelblue", "#FFBB00")
myColoursAlpha <- add.alpha(myColours, alpha=0.4)
my_cols <- setNames(myColoursAlpha, c("armature->no_armature", "no_armature->armature"))

pdf("../Fig4Pollination/armature_transition_times_mean_hist.pdf", height = 2.8, width = 6)
min <- min(avg_trans_times$avg_time)
max <- max(avg_trans_times$avg_time)
ax <- pretty(min:35, n = 20)

arm_to_no_arm <- avg_trans_times %>%
  dplyr::filter(transition == "armature->no_armature") %>%
  dplyr::mutate(arm_to_no_arm = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(arm_to_no_arm)

no_arm_to_arm <- avg_trans_times %>%
  dplyr::filter(transition == "no_armature->armature") %>%
  dplyr::mutate(no_arm_to_arm = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(no_arm_to_arm)

a2na <- hist(arm_to_no_arm$arm_to_no_arm, breaks = ax, plot = FALSE)
na2a <- hist(no_arm_to_arm$no_arm_to_arm, breaks = ax, plot = FALSE)

plot (a2na, col = myColours[1], xlab = "Time of transitions (mya)",  
      main = "", ylab = "number of transitions", 
      ylim = c(0, 40), xlim = c(100,0)) # alter if x values change!
plot (na2a, col = my_cols[2], add = TRUE)


# data for density graph
arm_to_no_arm <- armature_transitions %>%
  dplyr::filter(transition == "armature->no_armature") %>%
  dplyr::mutate(arm_to_no_arm = time) %>%
  dplyr::select(arm_to_no_arm, simulation = i)
no_arm_to_arm <- armature_transitions %>%
  dplyr::filter(transition == "no_armature->armature") %>%
  dplyr::mutate(no_arm_to_arm = time) %>%
  dplyr::select(no_arm_to_arm, simulation = i)

# calculate density curve
density_a2na <- density(arm_to_no_arm$arm_to_no_arm)
density_na2a <- density(no_arm_to_arm$no_arm_to_arm)

# plot the density
plot(density_a2na, lwd = 2, col = myColours[1], 
     xlim = c(58,0), xlab = "Time of transitions (mya)", bty = "l", ylim=c(0, 0.2),
     cex.lab = 1.4, cex.axis = 1.4, main = NULL, sub = NULL, title = NULL)
lines(density_na2a, lwd = 2, col = myColours[2], xlim = c(58,0))

# add data-points with noise in the X-axis
rug(jitter(no_arm_to_arm$no_arm_to_arm), col = my_cols[1])
rug(jitter(arm_to_no_arm$arm_to_no_arm), col = my_cols[2])

dev.off()
##################

pdf("Fig4Pollination/RatesPlot.pdf", height = 8, width = 10)
par(mar = c(5.1, 5.1, 4.1, 2.1))    # increase margins

arm_to_no_arm <- armature_trans_cumul %>%
  dplyr::filter(transition == "armature->no_armature")

# first set up basic plot parameters
plot(arm_to_no_arm$trans_no ~ arm_to_no_arm$time,
     type = "p", bty = "l", xlim = c(100,0), ylim = c(0,60),
     col = my_cols[1], pch = 15,
     xlab = "Time of transitions (mya)", 
     ylab = "Cumulative number of transitions",
     cex.lab = 1.8, cex.axis = 1.8)

# then loop through all data and add all points and lines to plot for 1000 simulations
for(n in 1:1000){
  test <- armature_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  no_arm_to_arm <- test %>%
    dplyr::filter(transition == "no_armature->armature")
  arm_to_no_arm <- test %>%
    dplyr::filter(transition == "armature->no_armature")
  
  points(arm_to_no_arm$trans_no ~ arm_to_no_arm$time, 
         col = my_cols[1], pch = 15)
  lines(arm_to_no_arm$trans_no ~ arm_to_no_arm$time, 
        col = my_cols[1])
  points(no_arm_to_arm$trans_no ~ no_arm_to_arm$time, 
         col = my_cols[2], pch = 17)
  lines(no_arm_to_arm$trans_no ~ no_arm_to_arm$time, 
        col =my_cols[2])
}

# then add average points and line in blue
# first prep data
no_arm_to_arm <- avg_trans_times %>%
  dplyr::filter(transition == "no_armature->armature") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
arm_to_no_arm <- avg_trans_times %>%
  dplyr::filter(transition == "armature->no_armature") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])

# then add to plot
points(arm_to_no_arm$trans_no ~ arm_to_no_arm$avg_time,
       col = "blue", pch = 15)
lines(arm_to_no_arm$trans_no ~ arm_to_no_arm$avg_time, 
      col = "blue")
points(no_arm_to_arm$trans_no ~ no_arm_to_arm$avg_time,
       col = "darkblue", pch = 17)
lines(no_arm_to_arm$trans_no ~ no_arm_to_arm$avg_time, 
      col = "darkblue")

# add legend
legend("topleft",
       legend = c("armature->no_armature", "no_armature->armature"),
       col = c(my_cols[1],my_cols[2]),
       pch = c(15, 17), pt.lwd = 0.001, bty = "n", cex = 1.8)
dev.off()
