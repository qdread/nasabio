# Fuzzed coordinates of trees

library(dplyr)

fia <- read.csv('/mnt/research/nasabio/data/fia/treedata10nov/finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv', stringsAsFactors = FALSE)

fiafuz <- fia %>%
	group_by(PLT_CN) %>%
	summarize(lat_fuzzed = na.omit(FUZZ_LAT)[1], lon_fuzzed = na.omit(FUZZ_LON)[1] )
	
library(sp)

aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fiasp <- SpatialPoints(coords = fiafuz[,c(3,2)], proj4string = CRS(wgs_crs))
fiaalbers <- spTransform(fiasp, CRSobj = CRS(aea_crs))

fiafuz <- cbind(fiafuz, lat_aea_fuzzed = fiaalbers@coords[,2], lon_aea_fuzzed = fiaalbers@coords[,1])
	
write.csv(fiafuz, '/mnt/research/nasabio/data/fia/fia_fuzzed_coords.csv', row.names = FALSE)