# Stratified random sampling of bbs or fia points, ensuring that the radii don't overlap, with a given radius size.
# QDR 03 Oct 2017

# Algorithm.
# Inputs: list of coordinates to subsample, polygon defining boundaries, radius around each point, and desired number of points.
# Make a settlers of catan hexagon grid on the landscape (to do this we need to select a starting point for offset and the appropriate size of each hexagon which will be determined by the buffer distance and possibly an additional parameter for the desired number of points
# Find a point in the interior of each hexagon that is not in the buffer perimeter zone. That's your subsample

# see strimas.com/spatial/hexagonal_grids

make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
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
                offset = c(0.5, 0.5))
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

library(maps)
usabounds <- map('usa')
