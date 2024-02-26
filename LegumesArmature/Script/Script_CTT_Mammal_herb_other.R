rm(list=ls())
library(geiger)
library(phytools)
library(dplyr)
library(evobiR)
library(RColorBrewer)


# Read in the data and tree ============================
##### list wd#####
trait <- ("C:\\Users\\rs67pisa\\Dropbox\\SEM_chapter1")
#####setwd mac####
trait01 <- ("/Users/rachelsouzaferreira/Dropbox/SEM_chapter1/Data")

setwd("C:\\Users\\rs67pisa\\Dropbox\\PhD\\Chapter 1\\SEMS_chapter1\\Data")
setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/SEMS_chapter1/Data")

#Load the trait data
trait <-read.csv("Trait_data.csv", sep=",")


##make dataframe with species, mass  trait
df1 <- trait[,c(1,6,12,20)]
subset_df <- df1[df1[, 2] == 1, ]

# Initialize a new column filled with zeros
subset_df$Condition_Met <- 0 

# Update this column to 1 for animals that meet both conditions (mass >= 10000 and diet >= 95)
subset_df$Condition_Met[subset_df$Mass.g >= 10000 & subset_df$Diet.Plant >= 95] <- 1


mammal_traits <- subset_df[,c(1,5)]
write.csv(mammal_traits, file = "mammal_traits.csv", row.names = FALSE)

str(mammal_traits)

library(geiger)
library(phytools)
library(ape)

#file.choose()

traitData <- read.csv("mammal_traits.csv", stringsAsFactors = FALSE, header = TRUE)
str(mammal_traits)
summary(mammal_traits)
dd<-mammal_traits 

setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Data/Phylacine/Phylogenies")
setwd("C:\\Users\\rs67pisa\\Dropbox\\PhD\\Chapter 1\\Data\\Phylacine\\Phylogenies")

# Load phylogeny
mammalPhylo <- read.nexus("Small_phylogeny.nex")
mammalPhylo2 <- read.nexus("Complete_phylogeny.nex")

tree <- mammalPhylo

# Load the required package
library(phangorn)

MCC <- mcc(tree, tree = TRUE, part = NULL, rooted = TRUE)

#save tree
write.nexus(tree, file="mammal_MCC.nex")



plot(MCC,cex = 0.5)


# Identify species that differ between the phylogeny and trait dataset
setdiff(MCC$tip.label, dd$Binomial.1.2)
setdiff(dd$Binomial.1.2, MCC$tip.label)

# Exclude those taxa from the phylogeny and the trait dataset
intersectTaxa <- intersect(dd$Binomial.1.2, MCC$tip.label)
traitDataSubset <- subset(dd, Binomial.1.2 %in% intersectTaxa)
mammalPhyloSubset <- drop.tip(MCC,
                              tip = MCC$tip.label[!MCC$tip.label %in% intersectTaxa])


# match data and tree ===================================
row.names(dd) <- dd$Binomial.1.2
setdiff(MCC$tip.label, dd$Binomial.1.2)


# Drop these species from the tree, because we need matching data between tree and traits
tree2 <- drop.tip(MCC, setdiff(MCC$tip.label, dd$Binomial.1.2))
tree2 <- ladderize(tree2)

matches <- match(dd$Binomial.1.2, tree2$tip.label, nomatch = 0)
data2 <- subset(dd, matches != 0)
data2 <- ReorderData(tree2, data2, taxa.names = "row names")
data2$Binomial.1.2 <- NULL
dd3 <- data2
str(dd3) # this data frame has the same species names as those in the tree and is ordered by the same order as the tip labels in the tree

# Select Mass trait (and make a vector from it)
Mass <- numeric(length = nrow(dd3))
unique(dd3$Condition_Met)
# Create the 'Mass' factor based on the 'Condition_Met' column
Mass <- ifelse(dd3$Condition_Met == 0, "other", "megaherbivore")
str(Mass)
head(Mass)
table(Mass)

# Convert to a factor
Mass <- as.factor(Mass)

# Display the levels of Mass to confirm
print(levels(Mass))

# Reorder data if needed (assuming ReorderData and tree2 are defined)
Mass <- ReorderData(tree2, Mass, taxa.names="names")
names(Mass) <- rownames(dd3)
# ordered trait vector for models (not dataframe)
str(Mass) # make sure it's a factor

tree <- tree2
trait <- Mass


setwd("C:\\Users\\rs67pisa\\Dropbox\\PhD\\Chapter 1\\Outputs")
setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Outputs")
# Clean environment to free memory:
rm(tree2,traitData, traitDataSubset, mammal_traits, dd)

# Mk models ============================
# estimate the rates of transitioning between ancestral states given a hypothesis (either equal rates(ER), or different rates(ARD))

# computing time = 2.297s
system.time(
  equal <- fitDiscrete(tree, trait, model = "ER"), gcFirst = T
)

# computing time = 8.726s
system.time(
  ard <- fitDiscrete(tree, trait, model = "ARD")
)

# Calculate AIC weights
aic.discrete <- setNames(c(equal$opt$aic, ard$opt$aic), c("equal", "different"))
weights <- aicw(aic.discrete)
weights # ARD model performs better for megaherbivore/no other


# For 100 simulation (nsim=100) --> computing time = 147.31s
system.time(
  mtrees <- make.simmap(tree, trait, model = "ARD", nsim = 100)
)

# 2.4s
system.time(
  pd_mtrees <- describe.simmap(mtrees, plot = FALSE)
)
# Save workspace in between =======================
save.image("workspace_after_Mk_models_MCC_mammals_herb.RData")
load(file = "Fig4Pollination_mammal/workspace_after_Mk_models_MCC.RData")

# We'd be interested in tracking the number of changes per unit of edge length rather than the total number of changes
# - because in all reconstructed phylogeny of extant taxa there is more edge length towards the tips of the tree than towards the root.
# we do this using the following argument: type = "rate"

### computing time = takes > 12h
#system.time(
 # obj <- ctt(mtrees[1], type="rate")
#)

###### Script from Pollination study: https://github.com/rubysaltbush/pollination-macroevolution/


# this function comes from https://github.com/rubysaltbush/pollination-macroevolution/blob/main/scripts/functions/transition_times.R

transition_times <- function(simmap){
  # below adapted from Liam Revells' phytools blog 
  # http://blog.phytools.org/2015/08/getting-timing-of-trait-changes-from.html
  # extracts raw transition times from a simmap (collapses multiple transitions
  # down into single transition events)
  # get tips and their states
  x <- phytools::getStates(simmap,"tips")
  # get unique states
  states <- sort(unique(x))
  # get length of states
  m <- length(states)
  # below makes a little matrix describing transitions
  ct <- sapply(states, 
               function(x,y) sapply(y, function(y,x) paste(x,"->", y, sep=""), 
                                    x = x), y = states)
  rm(x, states)
  # then a matrix to invalidate self->self transitions
  ii <- matrix(TRUE, m, m)
  diag(ii) <- rep(FALSE, m)
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
  # removing any nulls from list of changes and sorting other to large
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
mass_transitions <- data.frame()
for(i in 1:length(mtrees)){
  temp <- cbind(i, transition_times(mtrees[[i]]))
  mass_transitions <- rbind(mass_transitions, temp)
}
rm(temp, i)

table(mass_transitions$transition)

# build new data frame with cumulative number of transitions
mass_trans_cumul <- data.frame()
for(n in 1:1000){
  trans <- mass_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  mass_trans_cumul <- rbind(mass_trans_cumul, trans)
}
rm(n, trans)


##### copied from somewhere else in the code from the same script....:

# first need to know average number of transitions, rounded
trans_avg_length <- mass_trans_cumul %>%
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
avg_trans_times <- mass_trans_cumul %>%
  dplyr::group_by(transition, trans_no) %>%
  dplyr::mutate(avg_time = mean(time)) %>%
  dplyr::mutate(SE_time = sqrt(var(time) / length(time))) %>%
  dplyr::select(transition, trans_no, avg_time, SE_time) %>%
  dplyr::distinct()

# reduce avg_trans_times to trans_avg_length
avg_trans_times_a2na <- avg_trans_times %>%
  dplyr::filter(transition == "megaherbivore->other") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
avg_trans_times_na2a <- avg_trans_times %>%
  dplyr::filter(transition == "other->megaherbivore") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])
avg_trans_times <- rbind(avg_trans_times_a2na, avg_trans_times_na2a)
rm(avg_trans_times_a2na, avg_trans_times_na2a)

# export these results to csv in case I need them
readr::write_csv(avg_trans_times, "mean_transition_times_mass_MCC.csv")


#### Plotting ======================================

# install.packages("prettyGraphs") #to set transparancy of colors with "alpha"
library(prettyGraphs)
install.packages("prettyGraphs")

myColours = c("#FFBB00", "grey")
myColoursAlpha <- add.alpha(myColours, alpha=0.4)
my_cols <- setNames(myColoursAlpha, c("megaherbivore->other", "other->megaherbivore"))

pdf("mass_transition_times_mean_mammal_hist.pdf", height = 2.8, width = 6)
min <- min(avg_trans_times$avg_time)
max <- max(avg_trans_times$avg_time)

megaherbivore_to_other <- avg_trans_times %>%
  dplyr::filter(transition == "megaherbivore->other") %>%
  dplyr::mutate(megaherbivore_to_other = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(megaherbivore_to_other)

other_to_megaherbivore <- avg_trans_times %>%
  dplyr::filter(transition == "other->megaherbivore") %>%
  dplyr::mutate(other_to_megaherbivore = avg_time) %>%
  dplyr::ungroup() %>%
  dplyr::select(other_to_megaherbivore)

ax <- pretty(range(c(megaherbivore_to_other$megaherbivore_to_other, other_to_megaherbivore$other_to_megaherbivore)), n = 20)
a2na <- hist(megaherbivore_to_other$megaherbivore_to_other, breaks = ax, plot = FALSE)
na2a <- hist(other_to_megaherbivore$other_to_megaherbivore, breaks = ax, plot = FALSE)

# Adjust margins
par(mar = c(5, 5, 4, 2))

# Generate the plot with adjusted text size
plot(a2na, col = myColours[1], xlab = "Time of transitions (mya)",  
     main = "", ylab = "number of transitions", 
     ylim = c(0, 20), xlim = c(200, 0),
     cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8)

plot(na2a, col = my_cols[2], add = TRUE)

# data for density graph
megaherbivore_to_other <- mass_transitions %>%
  dplyr::filter(transition == "megaherbivore->other") %>%
  dplyr::mutate(megaherbivore_to_other = time) %>%
  dplyr::select(megaherbivore_to_other, simulation = i)
other_to_megaherbivore <- mass_transitions %>%
  dplyr::filter(transition == "other->megaherbivore") %>%
  dplyr::mutate(other_to_megaherbivore = time) %>%
  dplyr::select(other_to_megaherbivore, simulation = i)

# calculate density curve
density_a2na <- density(megaherbivore_to_other$megaherbivore_to_other)
density_na2a <- density(other_to_megaherbivore$other_to_megaherbivore)

# Calculate the range for the x-values
x_range_a2na <- range(density_a2na$x)
x_range_na2a <- range(density_na2a$x)

x_min <- min(x_range_a2na[1], x_range_na2a[1])
x_max <- max(x_range_a2na[2], x_range_na2a[2])

# Plot the density for a2na
plot(density_a2na, 
     lwd = 2, 
     col = myColours[1], 
     xlim = c(x_max, x_min), 
     xlab = "Time of transitions (mya)", 
     ylab = "Density", 
     main = NULL, 
     sub = NULL,
     bty = "l",
     cex.lab = 1.4, 
     cex.axis = 1.4)

# Add the density line for na2a
lines(density_na2a, 
      lwd = 2, 
      col = myColours[2])



# add data-points with noise in the X-axis
rug(jitter(other_to_megaherbivore$other_to_megaherbivore), col = my_cols[1])
rug(jitter(megaherbivore_to_other$megaherbivore_to_other), col = my_cols[2])

dev.off()
##################

pdf("RatesPlot_herbivore_other.pdf", height = 8, width = 10)
par(mar = c(5.1, 5.1, 4.1, 2.1))    # increase margins

megaherbivore_to_other<- mass_trans_cumul %>%
  dplyr::filter(transition == "megaherbivore->other")

# first set up basic plot parameters
plot(megaherbivore_to_other$trans_no ~ megaherbivore_to_other$time,
     type = "p", bty = "l", xlim = c(100,0), ylim = c(0,60),
     col = my_cols[1], pch = 15,
     xlab = "Time of transitions (mya)", 
     ylab = "Cumulative number of transitions",
     cex.lab = 1.8, cex.axis = 1.8)

# then loop through all data and add all points and lines to plot for 1000 simulations
for(n in 1:1000){
  test <- mass_transitions %>%
    dplyr::filter(i == n) %>%
    dplyr::group_by(transition) %>%
    dplyr::mutate(trans_no = row_number(i))
  
  other_to_megaherbivore<- test %>%
    dplyr::filter(transition == "other->megaherbivore")
  megaherbivore_to_other <- test %>%
    dplyr::filter(transition == "megaherbivore->other")
  
  points(megaherbivore_to_other$trans_no ~ megaherbivore_to_other$time, 
         col = my_cols[1], pch = 15)
  lines(megaherbivore_to_other$trans_no ~ megaherbivore_to_other$time, 
        col = my_cols[1])
  points(other_to_megaherbivore$trans_no ~ other_to_megaherbivore$time, 
         col = my_cols[2], pch = 17)
  lines(other_to_megaherbivore$trans_no ~ other_to_megaherbivore$time, 
        col =my_cols[2])
}

# then add average points and line in blue
# first prep data
other_to_megaherbivore <- avg_trans_times %>%
  dplyr::filter(transition == "other->megaherbivore") %>%
  dplyr::filter(trans_no <= trans_avg_length[1,2])
megaherbivore_to_other <- avg_trans_times %>%
  dplyr::filter(transition == "megaherbivore->other") %>%
  dplyr::filter(trans_no <= trans_avg_length[2,2])

# then add to plot
points(megaherbivore_to_other$trans_no ~ megaherbivore_to_other$avg_time,
       col = "blue", pch = 15)
lines(megaherbivore_to_other$trans_no ~ megaherbivore_to_other$avg_time, 
      col = "blue")
points(other_to_megaherbivore$trans_no ~ other_to_megaherbivore$avg_time,
       col = "darkblue", pch = 17)
lines(other_to_megaherbivore$trans_no ~ other_to_megaherbivore$avg_time, 
      col = "darkblue")

# add legend
legend("topleft",
       legend = c("other->megaherbivore", "megaherbivore->other"),
       col = c(my_cols[1],my_cols[2]),
       pch = c(15, 17), pt.lwd = 0.001, bty = "n", cex = 1.8)
dev.off()

