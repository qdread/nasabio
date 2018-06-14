# Stratified random sampling of bbs or fia points, ensuring that the radii don't overlap, with a given radius size.
# QDR 03 Oct 2017

# Algorithm.
# Inputs: list of coordinates to subsample, polygon defining boundaries, radius around each point, and desired number of points.
# Make a settlers of catan hexagon grid on the landscape (to do this we need to select a starting point for offset and the appropriate size of each hexagon which will be determined by the buffer distance and possibly an additional parameter for the desired number of points
# Find a point in the interior of each hexagon that is not in the buffer perimeter zone. That's your subsample

# see strimas.com/spatial/hexagonal_grids

make_hex_grid <- function(x, cell_diameter, cell_area, clip = FALSE, center_offset = c(0.5, 0.5)) {
  require(sp)
  require(raster)
  require(rgeos)
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate array of hexagon centers
  g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                offset = center_offset)
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
}

# For testing purposes, get a polygon with USA boundaries.
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

# Create hexagonal grid on the usa polygon, try 100km
usa_100k_hex <- make_hex_grid(x = usapoly, cell_diameter = 100000, clip = FALSE)

plot(usa_100k_hex)

# Set a buffer distance
# Then for each bbs centroid, determine which hexagon it is inside, and whether it falls in the buffer or not. (what is its distance from the polygon edge)
# If it falls in the buffer, there is a chance that its circle is overlapping other circles.

# Use a few random points as test points.
bbsll <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_correct_route_centroids.csv')

set.seed(8008135)
test_pts <- bbsll[sample(nrow(bbsll),10), ]

test_pts_spatial <- SpatialPoints(coords = test_pts[c('lon.1','lat.1')], proj4string = CRS(aea_crs))

inhex <- over(test_pts_spatial, usa_100k_hex) # This returns a vector saying which hexagon the point is in!

# For each point, determine the distance from the point to the center of the hexagon it's in.

# Get centroids
hex_centroids <- SpatialPoints(coords = getSpPPolygonsLabptSlots(usa_100k_hex), proj4string = CRS(aea_crs))

# Distance to centroids
centroid_dists <- sapply(1:10, function(i) gDistance(test_pts_spatial[i], hex_centroids[inhex[i]]))

# Convert polygons to lines
hex_asline <- as(usa_100k_hex, 'SpatialLines')

# Distance from point to edge of hexagon
edge_dists <- sapply(1:10, function(i) gDistance(test_pts_spatial[i], hex_asline[inhex[i]]))

# Any point that has an edge distance less than the buffer gets thrown out, then we see how many points are inside each hexagon.
# Let's use 25 km as the radius 
which(edge_dists >= 25000)
# 4 points are left

# SAMPLING FUNCTION
# Requires region and focal points to be in Albers projection (SpatialPolygons and SpatialPoints)
# Radius in m
# diameter in m
# Returns indices of the final points.
SRS_hexagons <- function(region, focal_points, radius, n, hex_diameter = NULL, offset = NULL) {
  
  # If diameter and offset aren't specified, create them.
  if (is.null(hex_diameter)) {
    hex_diameter <- radius * 2.75
  }
  if (is.null(offset)) {
    offset <- c(0.5, 0.5)
  }
  
  # Superimpose hexagonal grid on the region polygon
  hex_grid <- make_hex_grid(x = region, cell_diameter = hex_diameter, clip = FALSE, center_offset = offset)
  
  # Find which hexagon each point is in
  inhex <- over(focal_points, hex_grid)
  
  # Convert polygons to lines so that the distance is calculated properly
  hex_asline <- as(hex_grid, 'SpatialLines')
  
  # Distance from point to edge of hexagon
  # Add error checker because a few points might not be in a hexagon for some odd reason
  edge_dists <- sapply(1:length(focal_points), function(i) tryCatch(gDistance(focal_points[i], hex_asline[inhex[i]]), error = function(e) NA))
  
  # Throw out points that are too close to the edge of a hexagon
  use_points <- edge_dists >= radius & !is.na(edge_dists)
  
  # Subsample so that only one point per hexagon remains
  unique_hex <- unique(inhex[use_points])
  
  final_points <- sort(sapply(unique_hex, function(i) as.numeric(sample(names(inhex[inhex==i & use_points]), size=1))))
  
  # If there are less than n points, return the points with a warning.
  if (length(final_points) < n) {
    warning('Final number of points is less than n.')
    return(final_points)
  }
  
  # Otherwise, take a sample of size n of the remaining points.
  final_points <- sort(sample(final_points, size=n, replace=FALSE))
  return(final_points)
  
}


# Test function
bbs_aea <- SpatialPoints(coords=bbsll[,c('lon.1','lat.1')], proj4string = CRS(aea_crs))
SRS_hexagons(region=usapoly, focal_points = bbs_aea, radius = 25000, n = 1000, hex_diameter = 100000, offset = c(0.5, 0.5))
hexset1 <- SRS_hexagons(region=usapoly, focal_points = bbs_aea, radius = 25000, n = 500, hex_diameter = 100000, offset = c(0.5, 0.5))
hexset2 <- SRS_hexagons(region=usapoly, focal_points = bbs_aea, radius = 25000, n = 500, hex_diameter = 100000, offset = c(0.2, 0.7))

# See overlap between two sets
intersect(hexset1, hexset2) # Only 40 or so! Great!

# Try this with bigger ones though.
bighex1 <- SRS_hexagons(region=usapoly, focal_points=bbs_aea, radius=150000, n=500, offset=c(0.5,0.5))

# Try to do it with a lot of random offsets and different hexagon diameters

nsim <- 99
diams <- c(300e3, 400e3, 500e3, 600e3, 700e3)

sim_result <- matrix(NA, nrow=nsim, ncol=length(diams))

set.seed(8008)
pb <- txtProgressBar(0,nsim,style=3)

for (sim in 1:nsim) {
  setTxtProgressBar(pb,sim)
  for (diam in 1:length(diams)) {
    sim_result[sim, diam] <- length(SRS_hexagons(region=usapoly, focal_points=bbs_aea, radius=100e3, n=5000, hex_diameter=diams[diam], offset=runif(2)))
  }
}

close(pb)

# Clearly a small radius is better. Maybe even smaller?

nsim <- 10
diams <- c(250e3, 275e3)

sim_result <- matrix(NA, nrow=nsim, ncol=length(diams))

set.seed(8008008)
pb <- txtProgressBar(0,nsim,style=3)

for (sim in 1:nsim) {
  setTxtProgressBar(pb,sim)
  for (diam in 1:length(diams)) {
    sim_result[sim, diam] <- length(SRS_hexagons(region=usapoly, focal_points=bbs_aea, radius=100e3, n=5000, hex_diameter=diams[diam], offset=runif(2)))
  }
}

close(pb)