make_map_polygons <- function(states = character(0), countries = character(0)) {
  require(maps)
  require(maptools)
  require(raster)
  
  statepoly <- SpatialPolygons(list(), proj4string=CRS("+proj=longlat"))
  countrypoly <- SpatialPolygons(list(), proj4string=CRS("+proj=longlat"))
  
  if (length(states) > 0) {
    statebounds <- map('state', regions = tolower(states), fill = TRUE)
    stateIDs <- sapply(strsplit(statebounds$names, ":"), function(x) x[1])
    statepoly <- map2SpatialPolygons(statebounds, IDs=stateIDs, proj4string=CRS("+proj=longlat"))
    # Workaround to deal with the fact that Alaska isn't in the state database
    if ('alaska' %in% tolower(states)) {
      usabounds <- map('world', regions = 'usa', fill = TRUE)
      isak <- sapply(strsplit(usabounds$names, ":"), function(x) x[2])
      akpoly <- map2SpatialPolygons(usabounds, IDs=isak, proj4string=CRS("+proj=longlat"))['Alaska']
      statepoly <- bind(statepoly, akpoly)
    }
  }
  
  if (length(countries) > 0) {
    countrybounds <- map('world', regions = tolower(countries), fill = TRUE)
    countryIDs <- sapply(strsplit(countrybounds$names, ":"), function(x) x[1])
    countrypoly <- map2SpatialPolygons(countrybounds, IDs=countryIDs, proj4string=CRS("+proj=longlat"))
  }
  
  polys <- bind(statepoly, countrypoly)
  
}


# This function will flag coast plots
# Here coast plots are defined as anything near an ocean that is NOT near a land border
# So you need to first get the distance to the land borders for each point
# Then get the distance to the edge of the focal region for each point
# If the point is not close to the land border, but is close to the edge of the focal region, it's coast
flag_coast_plots <- function(focal_points, radius, focal_states = 'usa', border_states=character(0), border_countries=character(0)) {
  require(rgeos)
  aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
  # Create single polygon of focal region
  if (tolower(focal_states[1]) == 'usa') {
    focalpoly <- make_map_polygons(countries = 'usa')
  } else {
    focalpoly <- make_map_polygons(states=tolower(focal_states))
  }
  focalpoly <- gBuffer(focalpoly, width = 0)
  focalpoly <- as(focalpoly, 'SpatialLines')
  focalpoly <- spTransform(focalpoly, CRSobj = CRS(aea_crs))
  
  # Create polygons of bordering regions
  borderpoly <- make_map_polygons(states=border_states, countries=border_countries)
  borderpoly <- spTransform(borderpoly, CRSobj = CRS(aea_crs))
  
  # Test whether each of the sampling points has a radius that extends into any of the border areas
  # And test whether it's close to the edge of the focal area
  edge_distances <- c()
  coast_distances <- c()
  for (i in 1:length(focal_points)) {
    edge_distances[i] <- gDistance(focal_points[i], borderpoly)
    coast_distances[i] <- gDistance(focal_points[i], focalpoly)
  }
  
  # Return data frame with 3 columns, one is the distance, is_edge is whether it's an edge plot, and is_coast is whether it is a coast plot but not an edge plot.
  data.frame(edge_distance = edge_distances, 
             coast_distance = coast_distances, 
             is_edge = edge_distances < radius,
             is_coast = edge_distances >= radius & coast_distances < radius)
  
}
