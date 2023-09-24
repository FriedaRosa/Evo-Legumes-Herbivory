combined_ctt_plot <- function(trees, segments = 20, show.tree = FALSE, add = FALSE, ylim = NULL, xlim = NULL, lwd = 2, xlab = "time since the present", ylab = "mean change rate") {
  # Check if a tree belongs to the multiSimmap class
  if (!inherits(trees, "multiSimmap"))
    stop("trees should be an object of class \"multiSimmap\".")
  
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
  
  # Loop through the changes and segments
  for (i in 1:length(changes)) {
    for (j in 1:length(changes[[i]])) {
      for (k in 1:segments) {
        if ((changes[[i]][j] > segs[k, 1]) && (changes[[i]][j] <= segs[k, 2])) {
          nchanges[k] <- nchanges[k] + 1 / length(changes)
          edge.length[k] <- edge.length[k] + (min(segs[k, 2], changes[[i]][j]) - max(segs[k, 1], changes[[i]][j]))
        }
      }
    }
  }
  
  # Calculate the "rate" as changes per unit of edge length for each segment
  rate <- ifelse(edge.length > 0, nchanges / edge.length, 0)
  
  # Plot the "rate" type
  if (!add) {
    plot(segs[, 2], rate, type = "l", ylim = ylim, xlim = xlim, lwd = lwd, xlab = xlab, ylab = ylab)
  } else {
    lines(segs[, 2], rate, lwd = lwd)
  }
  
  # Optionally, add the tree to the plot
  if (show.tree) {
    plotTree(tree, add = TRUE, ftype = "off", lwd = 1, color = make.transparent("blue", 0.1), mar = par()$mar, direction = "leftwards", xlim = xlim)
  }
}

# Example usage:
combined_ctt_plot(trees, segments = 20, show.tree = TRUE)
