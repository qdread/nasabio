# Script to get HUC4 for each FIA plot.
# Author: QDR
# Project: NASA Bioxgeo
# Date: 14 Feb 2018

# Load true FIA plot coordinates.

fiacoords <- read.csv('~/data/allfia.csv')

# Load HUC4 polygon shape file.

fphuc <- '/mnt/research/aquaxterra/DATA/raw_data/HUC/shapefiles'

library(rgdal)
huc4 <- readOGR(dsn = fphuc, layer = 'WBDHU4')

# Run over() on the points and polygons.

fiacoords_sp <- SpatialPoints(coords = fiacoords[,c(3,2)],
							  proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
fiacoords_sp <- spTransform(fiacoords_sp,
							CRSobj = CRS(proj4string(huc4)))

fia_huc4 <- (fiacoords_sp %over% huc4)$HUC4 # Runs in < 5 minutes

# Combine output.
fia_huc <- data.frame(PLT_CN = fiacoords[,1],
                      HUC4 = fia_huc4)
					  
write.csv(fia_huc, '/mnt/research/nasabio/data/fia/fia_huc4.csv', row.names = FALSE)

#####

# Get already calculated HUC4 for each BBS route and make a CSV of it in the nasabio directory.

bbs_huc <- read.csv('/mnt/research/aquaxterra/CODE/python/BBSSpatialJoin/BBS_SpatialJoin_Final.csv', stringsAsFactors = FALSE)

# Find the one HUC4 containing the most stops for each BBS route and assign it there

library(dplyr)
library(tidyr)

bbs_huc4 <- bbs_huc %>%
	separate(rtestopNo, into = c('rteNo', 'Stop')) %>%
	group_by(rteNo) %>%
	summarize(HUC4 = names(sort(table(HUC4), decreasing = TRUE))[1])

write.csv(bbs_huc4, '/mnt/research/nasabio/data/bbs/bbs_huc4.csv', row.names = FALSE)