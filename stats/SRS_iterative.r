SRS_iterative <- function(focal_points, dist_mat = NULL, radius, n, point = NULL, show_progress = FALSE) {
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
    
    if(show_progress) print(paste('Found', s, 'of', n, 'points.'))
    
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

# Function that runs spDistsN1 with each iteration.
SRS_iterative_N1 <- function(focal_points, radius, n, point = NULL, show_progress = FALSE) {
  require(sp)
  
  if (is.null(point)) {
    point <- sample(length(focal_points), 1)
  }
  
  focal_points <- focal_points@coords
  dimnames(focal_points)[[1]] <- 1:nrow(focal_points)
  remaining_points <- focal_points
  
  final_points <- c()
  s <- 0
  
  # Keep going until the desired number of points is reached
  while(s < n) {
    final_points <- c(final_points, point)
    s <- s+1
    
    if(show_progress) print(paste('Found', s, 'of', n, 'points.'))
    
    dist_s <- spDistsN1(pts = remaining_points, pt = focal_points[dimnames(focal_points)[[1]] == point, ], longlat = FALSE)
    
    # Cross off all points within 2*radius from the focal point
    outside_circle <- dist_s >= 2*radius
    remaining_points <- remaining_points[outside_circle, , drop = FALSE]
    
    # If none are left, quit.
    if(nrow(remaining_points) == 0) {
      warning('Final number of points is less than n.')
      return(sort(as.numeric(final_points)))
    }
    
    # Find shortest distance in the remaining points to any one of the previously chosen points.
    point <- dimnames(remaining_points)[[1]][which.min(dist_s[outside_circle])]
  }
  
  return(sort(as.numeric(final_points)))
}
