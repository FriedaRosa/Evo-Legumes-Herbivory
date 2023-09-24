# Load the necessary packages
library(ape)
library(phytools)


## Function to get>Changes

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


# Simulate some changes on the tree
trees <- mtrees


# Define the number of segments
segments <- 20

# Check if a tree belongs to the multiSimmap class

# Convert the first tree to a phylo object
tree <- as.phylo(trees[[1]])

# Get the changes for all trees
changes <- sapply(trees, getChanges)

# Calculate the maximum height of the tree
h <- max(nodeHeights(tree))

# Create segment boundaries
segs <- cbind(seq(0, h - h / segments, h / segments),
              seq(1 / segments * h, h, h / segments))

# Initialize variables for calculating change rate
nchanges <- rep(0, segments)
edge.length <- rep(0, segments)











#### test

ii<-1
edge.length<-rep(0,b)
for(i in 1:nrow(segs)){
  done.seg<-FALSE
  while(LTT[ii,2]<=segs[i,2]&&done.seg==FALSE){
    print(i)
    edge.length[i]<-edge.length[i]+
      LTT[ii,1]*(min(segs[i,2],LTT[ii,3])-
                   max(segs[i,1],LTT[ii,2]))
    if(LTT[ii,3]>=segs[i,2]) done.seg<-TRUE
    if(LTT[ii,3]<=segs[i,2]) ii<-if(ii<nrow(LTT)) ii+1 else ii
  }
}

##########




# Loop through the changes and segments
for (i in 1:length(changes)) {                   #i = 1:100
  
  for (j in 1:length(changes[[i]])) {             #j = 1:137
    
    for (k in 1:segments) {                       #k = 1:20
      
      if ((changes[[i]][j] > segs[k, 1]) && (changes[[i]][j] <= segs[k, 2])) {
        
        nchanges[k] <- nchanges[k] + 1 / length(changes)
        edge.length[k] <- edge.length[k] + (min(segs[k, 2], changes[[i]][j]) - max(segs[k, 1], changes[[i]][j]))
        
      }
    }
  }
}

# Calculate the "rate" as changes per unit of edge length for each segment
rate <- nchanges / edge.length

# Create a plot of the "rate" type
plot(segs[, 2], rate, type = "l", ylim = c(0, max(rate)),
     xlim = c(max(segs[, 2]), min(segs[, 1])), lwd = 2,
     xlab = "time since the present", ylab = "mean change rate")
