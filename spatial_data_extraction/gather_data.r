# Put together all diversity values for Andy's analysis in wide format
# QDR NASAbioxgeo 14 Feb 2018

library(dplyr)

# FIA

# Geodiversity

fpdata <- '/mnt/research/nasabio/data'

# Point values
fiageopt <- read.csv(file.path(fpdata, 'fia/fia_geo_by_point.csv'), stringsAsFactors = FALSE)
fiahuc <- read.csv(file.path(fpdata, 'fia/fia_huc4.csv'), stringsAsFactors = FALSE)

fiadat <- fiageopt %>%
	select(PLT_CN, elevation_30m, bio1_1k, bio12_1k, geological_age_1k) %>%
	rename(elevation_point = elevation_30m, MAT_point = bio1_1k, MAP_point = bio12_1k, geoage_point = geological_age_1k) %>%
	left_join(fiahuc)
	
# Radius values

# Biodiversity

# Point values

# Radius values

###########################################################

# BBS

# Geodiversity

# Point values

bbsgeopt <- read.csv(file.path(fpdata, 'bbs/bbs_geo_by_point.csv'), stringsAsFactors = FALSE)
bbshuc <- read.csv(file.path(fpdata, 'bbs/bbs_huc4.csv'), stringsAsFactors = FALSE)

bbsdat <- bbsgeopt %>%
	select(rteNo, elevation_30m, bio1_1k, bio12_1k, geological_age_1k) %>%
	rename(elevation_point = elevation_30m, MAT_point = bio1_1k, MAP_point = bio12_1k, geoage_point = geological_age_1k) %>%
	left_join(bbshuc)

# Radius values
# Biodiversity
# Point values
# Radius values