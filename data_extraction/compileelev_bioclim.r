# Compile summary statistics: elevation and bioclim.
# BBS and FIA
# 04 Sep 2017 QDR

# Edited 13 Sep 2017: add geostats (footprint, soil type diversity, and geological age diversity)

library(dplyr)

##############################
# BBS

get_stats <- function(fp, prefix, n, suffix) {
	file_names <- paste0(prefix, n, suffix)
	all_stats <- list()
	for (i in file_names) {
		load(file.path(fp, i))
		all_stats[[length(all_stats) + 1]] <- stats_by_point
	}
	do.call('c', all_stats)
}

bbs_elev_stats <- get_stats('/mnt/research/nasabio/data/bbs/elevstats/30m', 'stats_', 1:1000, '.r')
bbs_bioclim5k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'bioclim5k_', 1:500, '.r')
bbs_bioclim1k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'bioclim1k_', 1:1000, '.r')
bbs_footprint_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'footprint_', 1:100, '.r')
bbs_geoage_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'geoage_', 1:100, '.r')
bbs_soil_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'soil_', 1:100, '.r')
bbs_night_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'nightlight_', 1:100, '.r')
bbs_dhi_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'dhi_', 1:100, '.r')

# Combine these, dealing with the fact that some list entries are NA.

bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv', stringsAsFactors = FALSE)
bbs_5k_missing <- bbs_bioclim5k_stats[[1]]
bbs_5k_missing[, 3:7] <- NA
bbs_1k_missing <- bbs_bioclim1k_stats[[1]]
bbs_1k_missing[, 3:7] <- NA
bbs_foot_missing <- bbs_footprint_stats[[1]]
bbs_foot_missing[, 3:7] <- NA
bbs_night_missing <- bbs_night_stats[[1]]
bbs_night_missing[, 3:7] <- NA
bbs_dhi_missing <- bbs_dhi_stats[[1]]
bbs_dhi_missing[, 3:7] <- NA


# geological age and soil type diversity have different column labels, so they will need to be joined differently.
bbs_geo_missing <- bbs_geoage_stats[[1]]
bbs_geo_missing[, 3:5] <- NA
bbs_soil_missing <- bbs_soil_stats[[1]]
bbs_soil_missing[, 3:5] <- NA

bbs_all_stats <- list()

for (i in 1:length(bbs_elev_stats)) {
	bbs_5k_i <- if (class(bbs_bioclim5k_stats[[i]]) == 'data.frame') bbs_bioclim5k_stats[[i]] else bbs_5k_missing
	bbs_5k_i$variable <- paste0('bio', 1:19, '_5k')
	bbs_1k_i <- if (class(bbs_bioclim1k_stats[[i]]) == 'data.frame') bbs_bioclim1k_stats[[i]] else bbs_1k_missing
	bbs_1k_i$variable <- paste0('bio', 1:19,'_1k')
	bbs_elev_i <- data.frame(radius = bbs_elev_stats[[i]][,'radius'], variable = 'elevation', bbs_elev_stats[[i]][,c('mean','sd','min','max')], n = NA)
	bbs_foot_i <- if(class(bbs_footprint_stats[[i]]) == 'data.frame') bbs_footprint_stats[[i]] else bbs_foot_missing
	bbs_foot_i$variable <- 'human_footprint'
	bbs_geo_i <- if(class(bbs_geoage_stats[[i]]) == 'data.frame') bbs_geoage_stats[[i]] else bbs_geo_missing
	bbs_geo_i$variable <- 'geological_age'
	bbs_soil_i <- if(class(bbs_soil_stats[[i]]) == 'data.frame') bbs_soil_stats[[i]] else bbs_soil_missing
	bbs_soil_i$variable <- 'soil_type'
	bbs_night_i <- if(class(bbs_night_stats[[i]]) == 'data.frame') bbs_night_stats[[i]] else bbs_night_missing
	bbs_night_i$variable <- 'nightlight'
	bbs_dhi_i <- if(class(bbs_dhi_stats[[i]]) == 'data.frame') bbs_dhi_stats[[i]] else bbs_dhi_missing
	bbs_dhi_i$variable <- rep(c('dhi_fpar', 'dhi_gpp', 'dhi_lai8', 'dhi_ndvi'), each = 8)
	
	bbs_all_stats[[i]] <- full_join(cbind(bbsll[i,], rbind(bbs_elev_i, bbs_5k_i, bbs_1k_i, bbs_foot_i, bbs_night_i, bbs_dhi_i)),
									cbind(bbsll[i,], rbind(bbs_geo_i, bbs_soil_i)))
}

bbs_all_stats <- do.call('rbind', bbs_all_stats)
bbs_all_stats <- rename(bbs_all_stats, lon_aea = lon.1, lat_aea = lat.1)

write.csv(bbs_all_stats, file = '/mnt/research/nasabio/data/bbs/bbs_geodiversity_stats.csv', row.names = FALSE)

##############################
# FIA

fia_elev_stats <- get_stats('/mnt/research/nasabio/data/fia/elevstats/newslice', 'stats_', 1:5000, '.r')
fia_bioclim5k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'bioclim5k_', 1:1000, '.r')
fia_bioclim1k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'bioclim1k_', 1:5000, '.r')
fia_footprint_stats <- get_stats('/mnt/research/nasabio/data/fia/geostats', 'footprint_', 1:250, '.r')
fia_geoage_stats <- get_stats('/mnt/research/nasabio/data/fia/geostats', 'geoage_', 1:250, '.r')
fia_soil_stats <- get_stats('/mnt/research/nasabio/data/fia/geostats', 'soil_', 1:250, '.r')
fia_night_stats <- get_stats('/mnt/research/nasabio/data/fia/geostats', 'nightlight_', 1:250, '.r')
fia_dhi_stats <- get_stats('/mnt/research/nasabio/data/fia/geostats', 'dhi_', 1:250, '.r')

# FIA lat long coordinates
library(dplyr)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

fia_5k_missing <- fia_bioclim5k_stats[[10001]]
fia_5k_missing[, 3:7] <- NA
fia_1k_missing <- fia_bioclim1k_stats[[10001]]
fia_1k_missing[, 3:7] <- NA
fia_elev_missing <- fia_elev_stats[[10001]]
fia_elev_missing[, 3:7] <- NA
fia_foot_missing <- fia_footprint_stats[[10001]]
fia_foot_missing[, 3:7] <- NA
fia_night_missing <- fia_night_stats[[10001]]
fia_night_missing[, 3:7] <- NA
fia_dhi_missing <- fia_dhi_stats[[10001]]
fia_dhi_missing[, 3:7] <- NA

# geological age and soil type diversity have different column labels, so they will need to be joined differently.
fia_geo_missing <- fia_geoage_stats[[10001]]
fia_geo_missing[, 3:5] <- NA
fia_soil_missing <- fia_soil_stats[[10001]]
fia_soil_missing[, 3:5] <- NA

fia_all_stats <- list()

for (i in 1:length(fia_elev_stats)) {
	fia_5k_i <- if (class(fia_bioclim5k_stats[[i]]) == 'data.frame') fia_bioclim5k_stats[[i]] else fia_5k_missing
	fia_5k_i$variable <- paste0('bio', 1:19)
	fia_1k_i <- if (class(fia_bioclim1k_stats[[i]]) == 'data.frame') fia_bioclim1k_stats[[i]] else fia_1k_missing
	fia_1k_i$variable <- paste0('bio', 1:19,'_1k')
	fia_elev_i <- if(class(fia_elev_stats[[i]]) == 'data.frame') fia_elev_stats[[i]] else fia_elev_missing
	fia_elev_i$variable <- 'elevation'
	fia_foot_i <- if(class(fia_footprint_stats[[i]]) == 'data.frame') fia_footprint_stats[[i]] else fia_foot_missing
	fia_foot_i$variable <- 'human_footprint'
	fia_geo_i <- if(class(fia_geoage_stats[[i]]) == 'data.frame') fia_geoage_stats[[i]] else fia_geo_missing
	fia_geo_i$variable <- 'geological_age'
	fia_soil_i <- if(class(fia_soil_stats[[i]]) == 'data.frame') fia_soil_stats[[i]] else fia_soil_missing
	fia_soil_i$variable <- 'soil_type'
	fia_night_i <- if(class(fia_night_stats[[i]]) == 'data.frame') fia_night_stats[[i]] else fia_night_missing
	fia_night_i$variable <- 'nightlight'
	fia_dhi_i <- if(class(fia_dhi_stats[[i]]) == 'data.frame') fia_dhi_stats[[i]] else fia_dhi_missing
	fia_dhi_i$variable <- rep(c('dhi_fpar', 'dhi_gpp', 'dhi_lai8', 'dhi_ndvi'), each = 13)
	

	fia_all_stats[[i]] <- full_join(cbind(fiall[i,], rbind(fia_elev_i, fia_5k_i, fia_1k_i, fia_foot_i, fia_night_i, fia_dhi_i)),
									cbind(fiall[i,], rbind(fia_geo_i, fia_soil_i)))
}			

fia_all_stats <- do.call('rbind', fia_all_stats)

write.csv(fia_all_stats, file = '/mnt/research/nasabio/data/fia/fia_geodiversity_stats.csv', row.names = FALSE)
