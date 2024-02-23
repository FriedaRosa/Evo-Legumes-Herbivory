#### Get Changes - 4 Traits Legumes ####

rm(list=ls())

####  Required functions ====================================================================== #
# ============================================================================================= #
# phytools-internal getChanges() function (part of ctt() function source code)
# ============================================================================================= #

getChanges<-function(tree){
  states<-sort(unique(getStates(tree)))
  nc<-sapply(tree$maps,length)-1
  b<-which(nc>0)
  nc<-nc[b]
  xx<-vector()
  H<-nodeHeights(tree)
  for(i in 1:length(b)){
    for(j in 1:nc[i]){
      ss<-names(tree$maps[[b[i]]])[j+1]
      x<-rep(H[b[i],1]+cumsum(tree$maps[[b[i]]])[j],2)
      xx<-c(xx,setNames(x[1],
                        paste(names(tree$maps[[b[i]]])[j:(j+1)],
                              collapse="->")))
    }
  }
  xx
}

# ============================================================================================= #
# Create a function to find the segment ID for a given time value:
# ============================================================================================= #

find_segment_id <- function(time) {
  segment_indices <- cut(time, breaks = segs_df[, 2], labels = FALSE)
  return(segs_df$segsID[segment_indices])
}

# ============================================================================ #


#### Data ===================================================================================== #
# ============================================================================================= #
# multiSimmap & pruned MCC tree matching the trait vector
# ============================================================================================= #


mtrees <- readRDS("../output/RDS/Legumes_4traits_simmap_100.rds") # multiSimmap
tree <- readRDS("../input/Mammals/Processed/Legumes_MCC_drop.rds") # MCC tree matching the data (pruned)


# Start:


# ================================================================= #
# Working the multiSiammp
# ================================================================= #

# Get all changes from all simulations so that we can loop through them:
changes <- sapply(mtrees, getChanges) 


# Transform the original changes list of matrices from sapply() for all simulations to dataframe 
# add trans_no to each transition based on grouping by transition 
# (i.e, there are two trans_no = 1, one for each transition type)


i_changes_list <- list() 
for (i in 1:length(changes)) {
  i_changes <- changes[[i]] # for each simulation from multiSimmap
  
  # extract transition types and times and assign ID column for simmap simulation
  i_changes_df <- data.frame(transition = names(i_changes), time = i_changes, N_sim = i)
  i_changes_df$transition <- as.factor(i_changes_df$transition)
  
  # Add trans_no to each transition from a group (transition type)
  i_changes_df <- i_changes_df %>%  arrange(time) %>%
    group_by(transition) %>%
    mutate(trans_no = row_number())
  
  i_changes_list[[i]] <- i_changes_df
}

# Final Table without edge.lengths:
changes_table <- plyr::rbind.fill(i_changes_list, fill=T)
changes_table_wo_edge <- changes_table

write.csv(changes_table_wo_edge, "../output/CSV/4traits_Legumes_changes_table_wo_edgelengths.csv")  
## This one works as well with the 4traits script, the edge-lengths are acrually not necessary, 
# but you can continue to run the code until the end if you want the edge.lengths to be part of the dataframe




# ====================================================================================================== #
# ~ Step 1: Loop over each simulation and extract 
# - breaks for segments
# - edge.lengths
# - rates for each character state transition
# ====================================================================================================== #

object_list <- list()
rates_list <- list()

# First Level of the loop  
for (Nsimulation in 1:100){ 
  
  # ====================================================================================================== #
  # Setting global variables for the Analysis
  # ====================================================================================================== #
  
  # set the number of time-segments for which everything will be calculated:
  b <- 20*5    
  
  # reduce to one tree to extract parameters from each simulation separately:
  tree <- as.phylo(mtrees[[Nsimulation]])
  
  # extract time since root
  h <- max(nodeHeights(tree))
  
  # calculate the time segments based on time since the origin and the time segments that were set by the user
  segs <- cbind(seq(0,h-h/b,h/b),
                seq(1/b*h,h,h/b))
  
  
  # ====================================================================================================== #
  # Transitions & Timing  (getChanges)
  # ====================================================================================================== #
  
  
  # ================================================================= #
  # Working changes from each simmap simulation separately:
  # ----------------------------------------------------------------- #
  # Changes per Segments:
  # ================================================================= #
  
  
  
  # Loop through all changes & Calculate average number of changes per time interval (segs; set with "b")
  
  nchanges <- rep(0,b)
  for(i in 1:length(changes)) { # loop through all simulations from simmap  
    
    # loop through all individual changes from one simmap simulation and calculate relative changes per time segment
    # nchanges = collects the relative number of changes (of all changes) per segment
    
    for(j in 1:length(changes[[i]])) { 
      ind<-which((changes[[i]][j]>segs[,1])+ 
                   (changes[[i]][j]<=segs[,2])==2)
      nchanges[ind]<-nchanges[ind]+1/length(changes)
      
    } # closing of 2nd level
  } # closing of 1st level
  
  rm(i,j)  
  
  # ================================================================= #
  # Working changes from each simmap simulation separately:
  # ----------------------------------------------------------------- #
  # Edge Lengths per Segments:
  # ================================================================= #
  
  # Edge.Lengths = Number of species in the phylogeny per time slot (increases naturally with time)
  # ! Note: has to be accounted for when showing transition rates
  # Here we calculate edge lengths for each segement separatly:
  
  # -------------------------------------------------------------------------------------------------------- #
  
  
  
  # Compute Lineages-Through-Time with ltt() function from the phytools package to extract edge.lengths from
  LTT <- ltt(tree,plot=FALSE, log.lineages = F)

  
  
  # transform LTT data so that we can use it:
  LTT <- cbind(LTT$ltt[2:(length(LTT$ltt)-1)],            # Number of species at each node minus the root
               LTT$times[2:(length(LTT$ltt)-1)],          # "start" time for each node minus the root
               LTT$times[3:length(LTT$ltt)])              # "end" time for each node minus the root
  
  
  ## Quick LTT plot coloured by trait
  
  # !!! Change trait states for armature here ========================================================= #
  # cols<-setNames(c("#009E50", "#709E01", "#E69F00", "#E68F90", "darkgrey"),
  #                c("herb_large","herb_small", "non_herb_large", "non_herb_small", "total")) 
  
  
  
  
  plot(LTT,col=cols,bty="n",las=1,lwd=4,show.tree=FALSE,ylim=c(1,Ntip(tree)),
       xlab="time (above the root)")
  # -------------------------------------------------------------------------------------------------#
  
  ## Now we loop through all segments to extract matching edge.lengths =================================== #
  
  # set count variable:
  ii<-1
  
  # Create an empty vector to store edge lengths
  edge.length <- rep(0, nrow(segs))
  
  
  # 2nd level loop: 
  # (for one tree because LTT is calculated for one tree)
  for (i_segs in 1:nrow(segs)) {
    
    
    # Find the indices of LTT segments that overlap with the current segs segment
    overlap_indices <- which(LTT[, 2] <= segs[i_segs, 2] & LTT[, 3] >= segs[i_segs, 1])
    
    
    # Calculate edge length for the current segs segment using vectorized operations
    # 3rd level loop: 
    for (ii in overlap_indices) {
      
      # edge.length should contain the edge lengths for all segments:  
      edge.length[i_segs] <- edge.length[i_segs] +        
        LTT[ii, 1] * # counts number of lineages between time intervals
        (min(segs[i_segs, 2], LTT[ii, 3]) - # end time interval
           max(segs[i_segs, 1], LTT[ii, 2]))  # start time interval
      
      
    } # 3rd level closing
    
  } # 2nd level closing
  
  
  edge.length_backup <- edge.length
  edge.length <-c(1, edge.length)
  # edge.length should contain the edge lengths for all segments in a given simulation
  
  # adding segsID coulmn to edge lengths:
  edge.length_df <- data.frame(edge.length = edge.length, segsID = seq(1,b+1))
  
  
  
  # ================================================================= #
  # Working the multiSimmap:
  # ----------------------------------------------------------------- #
  # Edge Lengths at each transition event
  # ================================================================= #
  ## Adapted code from the segments (above) by Frieda
  
  # we need the table for all simmaps: changes_table so that we can loop through it
  
  # set count variable:
  ii<-1
  
  # Create an empty vector to store edge lengths
  edge.length_transitions <- rep(0, nrow(changes_table))
  edge.lengths_all <- list()
  changes_times <- changes_table$time
  
  
  # 2nd level loop: 
  # (for one tree because LTT is calculated for one tree)
  for (i_changes in 1:nrow(changes_table)) {
    
    # Find the indices of LTT segments that overlap with the current segs segment
    
    # "which start times of the Lineages are smaller or equal to the time the change happened and which end times of the lineages are larger or equal to the time the change happend."
    # --> this gives the node numbers that were there when the transition happend
    
    overlap_indices <- which(LTT[, 2] <= changes_times[i_changes] & LTT[, 3] >= changes_times[i_changes])
    
    
    # Calculate edge length for the current segs segment using vectorized operations
    # 3rd level loop: 
    
    for (ii in overlap_indices) { # ii = number of lineages when a transition occurred
      
      # edge.length should contain the edge lengths for all segments:  
      edge.length_transitions[i_changes] <- edge.length_transitions[i_changes] +        
        LTT[ii, 1] 
      
      edge.length_transitions2 <- data.frame(time = changes_times[i_changes], edge.length = edge.length_transitions[i_changes])
      
      
      edge.lengths_all[[i_changes]] <- edge.length_transitions2
      
    } # 3rd level closing
    
  } # 2nd level closing
  
  
  edge.lengths_all_df <- plyr::rbind.fill(edge.lengths_all)
  changes_table_w_edge <- merge(changes_table_wo_edge, edge.lengths_all_df, by = "time")
  # Write to file (once is enough)
  write.csv(changes_table_w_edge, "../output/CSV/4traits_Legumes_changes_table_w_edgelengths.csv")