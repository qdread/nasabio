# Iterative search starting with a random point and finding the nearest point at least (2*radius) away from it that isn't already selected.
# Brute force method

# Calculate pairwise distance matrix for all the BBS points
library(sp)

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
bbsll <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_correct_route_centroids.csv')
bbs_aea <- SpatialPoints(coords=bbsll[,c('lon.1','lat.1')], proj4string = CRS(aea_crs))

bbs_dist <- spDists(bbs_aea, longlat = FALSE)

dimnames(bbs_dist) <- list(1:nrow(bbs_dist), 1:ncol(bbs_dist))

radius <- 100e3
n <- 200

# Start with random point
final_points <- c()
tosample <- rep(TRUE, length(bbs_aea))
point <- sample(length(bbs_aea), 1)
s <- 0

# Keep going until the desired number of points is reached
while(s < n) {
  final_points <- c(final_points, point)
  s <- s+1
  
  # Cross off all points within 2*radius from the focal point
  tosample[bbs_dist[point, ] < 2*radius] <- FALSE
  
  # If none are left, quit.
  if(!any(tosample)) break
  
  # Find shortest distance in the remaining points to any one of the previously chosen points.
  sub_dist <- bbs_dist[final_points, tosample, drop = FALSE]
  point <- as.integer(dimnames(sub_dist)[[2]][which(sub_dist == min(sub_dist), arr.ind = TRUE)[1,2]])
}

final_points <- sort(final_points)

# Plot them
library(sp)
library(maps)
library(maptools)
library(raster)
library(rgeos)

usabounds <- map('usa', fill=TRUE)
usaIDs <- sapply(strsplit(usabounds$names, ":"), function(x) x[1])
usapoly <- map2SpatialPolygons(usabounds, IDs = usaIDs, proj4string=CRS("+proj=longlat"))

# Transform USA polygon to a grid reference system, using our Albers projection.
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
usapoly <- spTransform(usapoly, CRSobj = CRS(aea_crs))

library(plotrix)

plot(usapoly)

for (i in 1:length(final_points)) {
  draw.circle(x = bbsll[final_points[i],'lon.1'], y = bbsll[final_points[i],'lat.1'], radius=100e3, border='skyblue', lwd=0.5)
}

points(bbs_aea[final_points], col = 'red', pch = 19, cex = 0.5)


#### SAMPLING FUNCTION
# Arguments: SpatialPoints object (and precalculated distance matrix if you want), radius, and sample size. The region does not matter. You can also supply a starting point index, but if you do not, a random one is provided.

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