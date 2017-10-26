SRS_iterative <- function(focal_points, dist_mat = NULL, radius, n, point = NULL) {
  require(sp)
  
  if (is.null(point)) {
    point <- sample(length(focal_points), 1)
  }
  if(is.null(dist_mat)) {
    dist_mat <- spDists(focal_points, longlat = FALSE)
  }
  
  dimnames(dist_mat) <- list(1:nrow(dist_mat), 1:ncol(dist_mat))
  
  final_points <- c()
  tosample <- rep(TRUE, length(focal_points))
  s <- 0
  
  # Keep going until the desired number of points is reached
  while(s < n) {
    final_points <- c(final_points, point)
    s <- s+1
    
    # Cross off all points within 2*radius from the focal point
    tosample[dist_mat[point, ] < 2*radius] <- FALSE
    
    # If none are left, quit.
    if(!any(tosample)) {
      warning('Final number of points is less than n.')
      return(sort(final_points))
    }
    
    # Find shortest distance in the remaining points to any one of the previously chosen points.
    sub_dist <- dist_mat[final_points, tosample, drop = FALSE]
    point <- as.integer(dimnames(sub_dist)[[2]][which(sub_dist == min(sub_dist), arr.ind = TRUE)[1,2]])
  }
  
  return(sort(final_points))
}