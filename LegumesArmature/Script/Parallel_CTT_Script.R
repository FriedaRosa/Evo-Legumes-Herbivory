# Load the necessary packages
library(doParallel)
library(foreach)

# Set the number of CPU cores to use
num_cores <- 4

# Set up a parallel backend with the specified number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

ctt_relative_parallel <- function(trees, segments = 20, num_cores = 4, ...) {
  if (!inherits(trees, "multiSimmap"))
    stop("trees should be an object of class \"multiSimmap\".")
  tree <- as.phylo(trees[[1]])
  h <- max(nodeHeights(tree))
  b <- segments
  segs <- cbind(seq(0, h - h / b, h / b), seq(1 / b * h, h, h / b))
  
  # Use foreach and %dopar% to parallelize the loop
  nchanges <- foreach(i = 1:length(trees), .combine = "c") %dopar% {
    changes <- sapply(trees[[i]], getChanges)
    nchanges_seg <- numeric(b)
    
    for (j in 1:length(changes)) {
      for (k in 1:length(changes[[j]])) {
        # Calculate relative changes
        ind <- which((changes[[j]][k] > segs[, 1]) + (changes[[j]][k] <= segs[, 2]) == 2)
        nchanges_seg[ind] <- nchanges_seg[ind] + 1 / length(changes)
      }
    }
    
    nchanges_seg
  }
  LTT <- ltt(tree, plot = FALSE)
  LTT <- cbind(LTT$ltt[2:(length(LTT$ltt) - 1)],
               LTT$times[2:(length(LTT$ltt) - 1)],
               LTT$times[3:length(LTT$ltt)])
  
  ii <- 1
  edge.length <- rep(0, b)
  
  for (i in 1:nrow(segs)) {
    done.seg <- FALSE
    while (LTT[ii, 2] <= segs[i, 2] && done.seg == FALSE) {
      edge.length[i] <- edge.length[i] +
        LTT[ii, 1] * (min(segs[i, 2], LTT[ii, 3]) -
                        max(segs[i, 1], LTT[ii, 2]))
      if (LTT[ii, 3] >= segs[i, 2]) done.seg <- TRUE
      if (LTT[ii, 3] <= segs[i, 2]) ii <- if (ii < nrow(LTT)) ii + 1 else ii
    }
  }
  
  # Return the result
  object <- list(segments = segs, nchanges = nchanges, edge.length = edge.length, tree = tree)
  class(object) <- "ctt_relative"
  object
  

}

# Close the parallel backend after using it
  stopCluster(cl)


tic()
# Load the necessary packages and define the ctt_relative_parallel function

# Set the number of CPU cores to use
num_cores <- 4

# Set up a parallel backend with the specified number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Use the ctt_relative_parallel function on your mtrees
result <- ctt_relative_parallel(mtrees, segments = 20)

# Close the parallel backend after using it
stopCluster(cl)

toc()