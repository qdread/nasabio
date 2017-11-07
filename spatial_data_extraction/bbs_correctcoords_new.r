# Read BBS shape files and calculate centroids for all routes.
# This will (1) make sure the most recent route data are included and (2) add back in the routes with loopy schwoopies.
# QDR 26 May 2017

library(rgdal)

# route data found at https://www.mbr-pwrc.usgs.gov/bbs/geographic_information/GIS_shapefiles_2015.html (downloaded 26 May 2017)
bbsrtes <- readOGR(dsn='C:/Users/Q/Dropbox/projects/nasabiodiv/bbsrtes_2012_alb', layer='bbsrte_2012_alb')

# Calculate centroids of each route.
# Route numbers should already be concatenated with state number.

# Update 31 May: need to combine the coordinates of routes that are represented by 2 segments into 1 segment.
library(dplyr)

centroid_all_segments <- function(x) {
  coord_x <- lapply(coordinates(x), '[[', 1)
  coord_x <- do.call(rbind, coord_x)
  cent_x <- apply(coord_x, 2, mean)
  data.frame(rteNo = x@data$rteno[1], lon = cent_x[1], lat = cent_x[2])
}

rteCentroids <- list()
for (rn in unique(bbsrtes@data$rteno)) {
  rteCentroids[[length(rteCentroids)+1]] <- centroid_all_segments(subset(bbsrtes, rteno == rn))
}

rteCentroids <- do.call(rbind, rteCentroids)

# Also transform to lat-long.
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'
rteCentroidsAlbers <- SpatialPoints(coords = rteCentroids[,2:3], proj4string = CRS(proj4string(bbsrtes)))
rteCentroidsLatLong <- spTransform(rteCentroidsAlbers, CRSobj = CRS(wgs_crs))

# Combine all route metadata and save as CSV.
bbs_coords <- cbind(rteNo = rteCentroids$rteNo, rteCentroidsLatLong@coords, rteCentroidsAlbers@coords)
names(bbs_coords) <- c('rteNo', 'lon','lat','lon_aea','lat_aea')
write.csv(bbs_coords, 'X:/data/bbs/bbs_correct_route_centroids.csv', row.names = FALSE)      
