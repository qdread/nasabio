# Read BBS shape files and calculate centroids for all routes.
# This will (1) make sure the most recent route data are included and (2) add back in the routes with loopy schwoopies.
# QDR 26 May 2017

library(rgdal)

# route data found at https://www.mbr-pwrc.usgs.gov/bbs/geographic_information/GIS_shapefiles_2015.html (downloaded 26 May 2017)
# alternatively see https://earthworks.stanford.edu/catalog/stanford-vy474dv5024 (downloaded 17 Apr 2018)
bbsrtes <- readOGR(dsn='C:/Users/Q/Dropbox/projects/nasabiodiv/bbsrtes_2012_alb', layer='bbsrte_2012_alb')
bbsrtes2014 <- readOGR(dsn='C:/Users/Q/Dropbox/projects/nasabiodiv/bbsrtes2014', layer='bbsrtsl020')

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

    

# Find route midpoints ----------------------------------------------------

# Step One, find routes that only consist of a single segment
# Step Two, find the midpoint of each of those segments

# 4111 routes, 4740 segments

library(dplyr)
library(maptools)

n_segments <- table(bbsrtes@data$rteno)
table(n_segments) # 3754 are in a single piece

single_segments <- bbsrtes@data %>%
  group_by(rteno) %>%
  summarize(n_segments = n()) %>%
  filter(n_segments == 1)

bbsrtes_single <- subset(bbsrtes, rteno %in% single_segments$rteno)
bbsrtes_mid <- SpatialLinesMidPoints(bbsrtes_single)

# Three routes still have multiple midpoints after doing this, so get rid of them to go down to 3751 routes
bad_routes <- bbsrtes_mid@data %>%
  group_by(rteno) %>%
  summarize(n_segments = n()) %>%
  filter(n_segments > 1)

bbsrtes_mid <- subset(bbsrtes_mid, !rteno %in% bad_routes$rteno)
bbsrtes_mid_latlong <- spTransform(bbsrtes_mid, CRSobj = CRS('+proj=longlat +ellps=WGS84 +no_defs'))

bbsrtes_mid_df <- data.frame(rteNo = bbsrtes_single@data$rteno)
bbsrtes_mid_lldf <- data.frame(rteNo = bbsrtes_mid@data$rteno,
                               lon = bbsrtes_mid_latlong@coords[,1],
                               lat = bbsrtes_mid_latlong@coords[,2],
                               lon_aea = bbsrtes_mid@coords[,1],
                               lat_aea = bbsrtes_mid@coords[,2])

bbsrtes_mid_df <- left_join(bbsrtes_mid_df, bbsrtes_mid_lldf) %>% filter(complete.cases(.))
write.csv(bbsrtes_mid_df, 'C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_route_midpoints.csv', row.names = FALSE)


# Plot locations of routes with multiple segments
broken_bbsrtes <- subset(bbsrtes, !rteno %in% bbsrtes_mid_df$rteNo)
broken_bbsrtes <- spTransform(broken_bbsrtes, CRSobj = CRS('+proj=longlat +ellps=WGS84 +no_defs'))

library(maps)
map('state')
lines(broken_bbsrtes, col = 'red', main = 'Locations of non-continuous BBS routes')

# There are a lot of them in western NC.
map('state', xlim = c(-84.4, -80), ylim = c(34, 37))
lines(broken_bbsrtes, col = 'red', main = 'Locations of non-continuous BBS routes')

# Save table of bbs routes with multiple segments
broken_bbsrtes@data %>%
  select(rteno, RTENAME) %>%
  unique %>%
  write.csv('C:/Users/Q/google_drive/NASABiodiversityWG/DataPaper/noncontinuousbbsroutes.csv', row.names = FALSE)
