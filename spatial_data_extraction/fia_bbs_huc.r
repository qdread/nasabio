# Script to get HUC4 for each FIA plot.
# Author: QDR
# Project: NASA Bioxgeo
# Date: 14 Feb 2018

# Edited 25 April 2018: use BBS midpoints, not centroids.
# Edited 05 March 2018: add Bird Conservation Regions (BCR) and TNC Ecoregions (TNC)

# Load true FIA plot coordinates.

fiacoords <- read.csv('~/data/allfia.csv')

# Load HUC4 polygon shape file.

fphuc <- '/mnt/research/aquaxterra/DATA/raw_data/HUC/shapefiles'
fpregion <- '/mnt/research/nasabio/data/ecoregions'

library(rgdal)
huc4 <- readOGR(dsn = fphuc, layer = 'WBDHU4')
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
# TNC must be converted to WGS84
tnc <- spTransform(tnc, CRSobj = CRS(proj4string(huc4)))

# Run over() on the points and polygons.

fiacoords_sp <- SpatialPoints(coords = fiacoords[,c(3,2)],
							  proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
fiacoords_sp <- spTransform(fiacoords_sp,
							CRSobj = CRS(proj4string(huc4)))

fia_huc4 <- (fiacoords_sp %over% huc4)$HUC4 # Runs in < 5 minutes
fia_bcr <- (fiacoords_sp %over% bcr)$BCRNAME
fia_tnc <- (fiacoords_sp %over% tnc)$ECODE_NAME

# Combine output.
fia_huc <- data.frame(PLT_CN = fiacoords[,1],
                      HUC4 = fia_huc4,
					  BCR = fia_bcr,
					  TNC = fia_tnc)
					  
write.csv(fia_huc, '/mnt/research/nasabio/data/fia/fia_ecoregions.csv', row.names = FALSE)

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

# Load lat long coordinates for centroid of each route and find the BCR and TNC it's in

#bbscoords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv', stringsAsFactors = FALSE)
bbscoords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv', stringsAsFactors = FALSE)
bbscoords_sp <- SpatialPoints(coords = bbscoords[,c('lon', 'lat')],
							  proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
bbscoords_sp <- spTransform(bbscoords_sp,
							CRSobj = CRS(proj4string(huc4)))

bbs_bcr <- (bbscoords_sp %over% bcr)$BCRNAME
bbs_tnc <- (bbscoords_sp %over% tnc)$ECODE_NAME

bbs_ecoregions <- read.csv('/mnt/research/nasabio/data/bbs/bbs_huc4.csv', stringsAsFactors = FALSE) %>%
	left_join(data.frame(rteNo = bbscoords$rteNo,
						 BCR = bbs_bcr,
						 TNC = bbs_tnc))
						 
write.csv(bbs_ecoregions, '/mnt/research/nasabio/data/bbs/bbs_ecoregions.csv', row.names = FALSE)

# Load bbs all geodiversity data frame and join the ecoregions to it (in right order)
### this section is no longer needed 25 apr 2018
# bbsgeo <- read.csv('/mnt/research/nasabio/data/bbs/bbs_allgeo_wide.csv', stringsAsFactors = FALSE)
# bbsgeo <- bbsgeo %>%
	# left_join(bbs_ecoregions) %>%
	# select(rteNo, HUC4, BCR, TNC, everything())
	
# write.csv(bbsgeo, '/mnt/research/nasabio/data/bbs/bbs_allgeo_wide.csv', row.names = FALSE)