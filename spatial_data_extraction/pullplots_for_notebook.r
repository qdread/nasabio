# Identify which BBS route midpoints, and FIA plot locations, are within 100 km of the state borders of New Hampshire.
# QDR/NASAbioXgeo/09 Oct 2018

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0' # for NLCD
wgs_crs <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' # for elevation layers


library(dplyr)
library(maps)
library(maptools)
library(rgeos)
library(data.table)

bbscoords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv', stringsAsFactors = FALSE) %>% select(rteNo, lat, lon)
fiadat <- read.csv('/mnt/research/nasabio/data/fia/treedata10nov/finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv', 
				   stringsAsFactors = FALSE, 
				   colClasses = c(STATECD = 'character', PLT_CN = 'character', TREE_CN = 'character', UNITCD = 'character', COUNTYCD = 'character'))

# 100 km buffer around New Hampshire.				   
nhpoly_albers_buffer <- 
	map('state', regions = 'new hampshire', fill = TRUE) %>%
	map2SpatialPolygons(IDs = 'new hampshire', proj4string=CRS("+proj=longlat")) %>%
	spTransform(CRSobj = CRS(aea_crs)) %>%
	gBuffer(width = 100e3)

# Identify which BBS midpoints fall within the buffer.	
bbs_inbuffer <- bbscoords %>%
	select(lon, lat) %>%
	SpatialPoints(proj4string = CRS(wgs_crs)) %>%
	spTransform(CRSobj = CRS(aea_crs)) %>%
	over(nhpoly_albers_buffer) 
	
bbsrtes_inbuffer <- bbscoords$rteNo[!is.na(bbs_inbuffer)]

# FIPS codes for New England/NY (quick subset)
ne_ny <- c(9, 23, 25, 33, 36, 44, 50)

fia_inbuffer <- fiadat %>%
	filter(STATECD %in% ne_ny) %>%
	select(FUZZ_LON, FUZZ_LAT) %>%
	SpatialPoints(proj4string = CRS(wgs_crs)) %>%
	spTransform(CRSobj = CRS(aea_crs)) %>%
	over(nhpoly_albers_buffer) 
	
fiaplots_inbuffer <- fiadat %>%
	filter(STATECD %in% ne_ny) %>%
	filter(!is.na(fia_inbuffer)) %>%
	pull(PLT_CN) %>%
	unique

# Create subset of BBS data from "raw" data
fpbbs <- '/mnt/research/aquaxterra/DATA/raw_data/BBS/DataFiles26may2017/50-StopData/1997ToPresent_SurveyWide'
bbsdat <- lapply(1:10, function(i) read.csv(file.path(fpbbs, paste0('fifty', i, '.csv')), stringsAsFactors = FALSE, colClasses = c(statenum = 'character', Route = 'character')))
bbsdat <- do.call(rbind, bbsdat)
bbsdat <- bbsdat %>%
	mutate(rteNo = paste(statenum, Route, sep = ''))

bbs_subset <- bbsdat %>%
	filter(rteNo %in% bbsrtes_inbuffer, between(year, 2007, 2016))

bbs_subset <- bbs_subset %>%
	select(-rteNo)
	
# Create subset of FIA data from "raw" data.
# Also pre-remove the ones that are not plantation plots.
plotcond <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv', colClasses = c(PLT_CN = 'character'))

fia_subset <- fiadat %>%
	filter(PLT_CN %in% fiaplots_inbuffer, PLT_CN %in% plotcond$PLT_CN[!plotcond$plantation])
	
# Write both the CSVs.
write.csv(bbs_subset, '/mnt/research/nasabio/temp/sampledata/bbs_data.csv', row.names = FALSE)	
write.csv(fia_subset, '/mnt/research/nasabio/temp/sampledata/fia_data.csv', row.names = FALSE)	