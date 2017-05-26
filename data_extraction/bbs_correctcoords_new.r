# Read BBS shape files and calculate centroids for all routes.
# This will (1) make sure the most recent route data are included and (2) add back in the routes with loopy schwoopies.
# QDR 26 May 2017

library(rgdal)

# route data found at https://www.mbr-pwrc.usgs.gov/bbs/geographic_information/GIS_shapefiles_2015.html (downloaded 26 May 2017)
bbsrtes <- readOGR(dsn='C:/Users/Q/Dropbox/projects/nasabiodiv/bbsrtes_2012_alb', layer='bbsrte_2012_alb')

# Calculate centroids of each route.
# Route numbers should already be concatenated with state number.
rteCentroids <- lapply(coordinates(bbsrtes), function(x) apply(x[[1]], 2, mean))
rteCentroids <- do.call(rbind, rteCentroids)

# Also transform to lat-long.
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'
rteCentroidsAlbers <- SpatialPoints(coords = rteCentroids, proj4string = CRS(proj4string(bbsrtes)))
rteCentroidsLatLong <- spTransform(rteCentroidsAlbers, CRSobj = CRS(wgs_crs))

# Combine all route metadata and save as CSV.
bbs_coords <- cbind(bbsrtes@data, rteCentroidsLatLong@coords, rteCentroidsAlbers@coords)
names(bbs_coords)[7:10] <- c('lon','lat','lon_aea','lat_aea')
names(bbs_coords)[1] <- 'rteNo'
write.csv(bbs_coords, 'X:/data/bbs/bbs_correct_route_centroids.csv', row.names = FALSE)      
