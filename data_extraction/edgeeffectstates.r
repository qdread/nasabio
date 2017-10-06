# Get rid of nasabioxgeo plots (bbs or fia) that have radii that extend into neighboring land areas that don't have biodiversity measured
# This will deal with some of the beta-diversity edge effects
# For the Pacific Northwest analysis, we will have to get rid of anything along the land border of CA, OR, and WA with other states and with Can/Mex
# For the entire USA analysis, we will only have to get rid of ones that extend into Can/Mex
# It's possible to get rid of coasts too, but that's a different kind of edge effect.


# Procedure:
# 1. Get polygons of the actual sampling area (PNW or USA) and polygons of the non-sampled area (Can/Mex/ID/NV/AZ/UT/MT; Can/Mex only)
# 2. Define radius r
# 3. Find which points are within r distance of the non-sampled area and flag them.

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

flag_edge_plots <- function(focal_points, radius, border_states=character(0), border_countries=character(0)) {
  require(rgeos)
  
  # Create polygons of bordering regions
  borderpoly <- make_map_polygons(states=border_states, countries=border_countries)
  aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
  borderpoly <- spTransform(borderpoly, CRSobj = CRS(aea_crs))
  
  # Test whether each of the sampling points has a radius that extends into any of the polygons
  distances <- c()
  for (i in 1:length(focal_points)) {
    distances[i] <- gDistance(focal_points[i], borderpoly)
  }
  
  # Return data frame with 2 columns, one is the distance, and the other is whether it's an edge plot
  data.frame(distance = distances, is_edge = distances <= radius)
  
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
  
  # Return data frame with 2 columns, one is the distance, and the other is whether it's an edge plot
  data.frame(edge_distance = edge_distances, 
             coast_distance = coast_distances, 
             is_edge = edge_distances < radius,
             is_coast = edge_distances >= radius & coast_distances < radius)
  
}

# Transform USA polygon to a grid reference system, using our Albers projection.
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'


states <- c('Arizona','Nevada','Idaho')
countries <- c('Canada','Mexico')

library(dplyr)
fiapnw <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
#fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

fia_aea <- spTransform(SpatialPoints(coords=fiacoords[,c('lon','lat')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))
gDistance(fia_aea,borderpoly)

fia_edge100 <- flag_edge_plots(focal_points = fia_aea, radius = 100e3, border_states = c('Arizona','Idaho','Nevada'), border_countries = c('Canada', 'Mexico'))

borderpoly_fia <- make_map_polygons(states=c('Arizona','Idaho','Nevada'), countries=c('Canada', 'Mexico'))
borderpoly_fia <- spTransform(borderpoly_fia, CRSobj = CRS(aea_crs))

plot(borderpoly_fia)
points(fia_aea, col = c('red','blue')[as.numeric(fia_edge100$is_edge) + 1], pch = 19, cex = 0.3)


bbsll <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_correct_route_centroids.csv')
bbs_aea <- SpatialPoints(coords=bbsll[,c('lon.1','lat.1')], proj4string = CRS(aea_crs))

bbs_edge100 <- flag_edge_plots(focal_points = bbs_aea, radius = 100e3, border_countries = c('Canada', 'Mexico'))

borderpoly_bbs <- make_map_polygons(countries = c('Canada','Mexico'))
borderpoly_bbs <- spTransform(borderpoly_bbs, CRSobj = CRS(aea_crs))

plot(borderpoly_bbs)
points(bbs_aea, col = c('red','blue')[as.numeric(bbs_edge100$is_edge) + 1], pch=19, cex=0.3)

fia_coast100 <- flag_coast_plots(focal_points = fia_aea, radius = 100e3, focal_states = c('California','Oregon','Washington','Alaska'), border_states = c('Arizona','Idaho','Nevada'), border_countries = c('Canada', 'Mexico'))

focalpoly_fia <- make_map_polygons(states = c('California','Oregon','Washington'))
focalpoly_fia <- gBuffer(focalpoly_fia, width = 0)
focalpoly_fia <- spTransform(focalpoly_fia, CRSobj = CRS(aea_crs))

plot(borderpoly_fia)
points(fia_aea, col = c('red','blue')[as.numeric(fia_coast100$is_coast) + 1], pch = 19, cex = 0.3)

plot(focalpoly_fia)
points(fia_aea[!fia_coast100$is_edge & !fia_coast100$is_coast], pch=19, cex=0.3, col='green')
points(fia_aea[fia_coast100$is_edge & !fia_coast100$is_coast], pch=19, cex=0.3, col='red')
points(fia_aea[!fia_coast100$is_edge & fia_coast100$is_coast], pch=19, cex=0.3, col='blue')
lines(gBuffer(focalpoly_fia, width = -100e3))