# compare bbs coordinates in my spreadsheet to the "real" ones given by BBS data file.

fp <- '/mnt/research/nasabio/data/bbs'
routes <- read.csv(file.path(fp,'routes.csv'), stringsAsFactors=FALSE)
oldcoords <- read.csv(file.path(fp, 'bbs_wgs84_coords_byroute.csv'), stringsAsFactors=FALSE)

# paste together route names.
# add leading zeroes.
firstzero <- c('','0')[as.numeric(routes$Route < 100) + 1]
secondzero <- c('','0')[as.numeric(routes$Route < 10) + 1]

routes$rteNo <- as.numeric(paste0(as.numeric(routes$statenum), firstzero, secondzero, routes$Route))
write.csv(routes, file.path(fp, 'bbs_correct_coords_byroute.csv'), row.names=FALSE)


######################################
# unfortunately the ones above are the starting point. We want the centroid.

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
latlong_crs <- '+proj=longlat +ellps=WGS72 +no_defs'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fp <- '/mnt/research/nasabio/data/bbs'
bbsll <- read.csv(file.path(fp, 'bbs_wgs72_correct_coords.csv'), stringsAsFactors = FALSE)

# Get route centroids in latitude and longitude (ignoring year)
# Convert them to Albers.

library(dplyr)

routeCentroid <- function(x) {
	stopcoords <<- unique(cbind(x$POINT_X, x$POINT_Y))
	data.frame(lon = mean(stopcoords[,1]), lat = mean(stopcoords[,2]))
}

bbs_centroids <- bbsll %>% 
	filter(!is.na(POINT_X)) %>%
	group_by(rteNo) %>% 
	do(routeCentroid(.))
	
library(sp)
bbs_centroids <- as.data.frame(bbs_centroids)
bbs_centroids_albers <- spTransform(SpatialPointsDataFrame(coords=bbs_centroids[,2:3], data=bbs_centroids[,1,drop=FALSE], proj4string=CRS(latlong_crs)), CRSobj=CRS(aea_crs))
bbs_centroids_wgs84 <- spTransform(SpatialPointsDataFrame(coords=bbs_centroids[,2:3], data=bbs_centroids[,1,drop=FALSE], proj4string=CRS(latlong_crs)), CRSobj=CRS(wgs_crs))
# The lat longs with different reference ellipsoids are so similar that it does not matter for our purposes which one is used.

bbs_centroids <- cbind(bbs_centroids, lon_aea = bbs_centroids_albers@coords[,1], lat_aea = bbs_centroids_albers@coords[,2])
write.csv(bbs_centroids, file.path(fp, 'bbs_correct_route_centroids.csv'), row.names = FALSE)