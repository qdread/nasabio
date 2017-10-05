# Compile summary statistics: elevation and bioclim.
# BBS and FIA
# 04 Sep 2017 QDR

# Edited 04 Oct 2017: add new elevation stats that include elevation, slope, sin aspect, cos aspect, and TPI
# Edited 25 Sep 2017: add biocloud
# Edited 13 Sep 2017: add geostats (footprint, soil type diversity, and geological age diversity)

library(dplyr)

##############################
# BBS

get_stats <- function(fp, prefix, n, suffix, stacked = FALSE) {
	file_names <- paste0(prefix, n, suffix)
	all_stats <- list()
	for (i in file_names) {
		load(file.path(fp, i))
		if (!stacked) all_stats[[length(all_stats) + 1]] <- stats_by_point
		if (stacked) all_stats[[length(all_stats) + 1]] <- lapply(stats_by_point, function(x) do.call('rbind', x))
	}
	do.call('c', all_stats)
}

get_all_elev_stats <- function(fp, prefix, n, suffix) {
	file_names <- paste0(prefix, n, suffix)
	elev_stats <- list()
	sin_aspect_stats <- list()
	cos_aspect_stats <- list()
	slope_stats <- list()
	TPI_stats <- list()
	for (i in file_names) {
		load(file.path(fp, i))
		elev_stats[[length(elev_stats) + 1]] <- lapply(stats_by_point, '[[', 1)
		slope_stats[[length(slope_stats) + 1]] <- lapply(stats_by_point, '[[', 3)
		TPI_stats[[length(TPI_stats) + 1]] <- lapply(stats_by_point, '[[', 4)
		sin_aspect_stats[[length(sin_aspect_stats) + 1]] <- tryCatch(lapply(stats_by_point, function(x) {
			dat <- x[[2]][,c('radius','variable','mean_sin','sd_sin','min_sin','max_sin','n')]
			names(dat) <- c('radius','variable','mean','sd','min','max','n')
			dat
		}), error = function(e) rep(NA, length(stats_by_point)))
		cos_aspect_stats[[length(cos_aspect_stats) + 1]] <- tryCatch(lapply(stats_by_point, function(x) {
			dat <- x[[2]][,c('radius','variable','mean_cos','sd_cos','min_cos','max_cos','n')]
			names(dat) <- c('radius','variable','mean','sd','min','max','n')
			dat
		}), error = function(e) rep(NA, length(stats_by_point)))
	}
	list(elev = do.call('c', elev_stats),
		 slope = do.call('c', slope_stats),
		 TPI = do.call('c', TPI_stats),
		 sin_aspect = do.call('c', sin_aspect_stats),
		 cos_aspect = do.call('c', cos_aspect_stats))
}

replace_na_df <- function(dflist) {
	good <- sapply(dflist, is.data.frame)
	missing_df <- dflist[[which(good)[1]]]
	missing_df[, 3:ncol(missing_df)] <- NA
	lapply(dflist, function(x) if (!inherits(x, 'data.frame')) missing_df else x)
}

replace_varname <- function(dflist, varname) {
	lapply(dflist, function(x) {
		x$variable <- varname
		x
	})
}

bbs_elev_stats <- get_all_elev_stats('/mnt/research/nasabio/data/bbs/elevstats/big30m', 'stats_', 1:2000, '.r')
bbs_elevation_stats <- bbs_elev_stats$elev
bbs_slope_stats <- bbs_elev_stats$slope
bbs_tpi_stats <- bbs_elev_stats$TPI
bbs_sin_aspect_stats <- bbs_elev_stats$sin_aspect
bbs_cos_aspect_stats <- bbs_elev_stats$cos_aspect

bbs_bioclim5k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'bioclim5k_', 1:500, '.r')
bbs_bioclim1k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'bioclim1k_', 1:1000, '.r')
bbs_biocloud5k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'biocloud5k_', 1:1000, '.r', stacked = TRUE)
bbs_biocloud1k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'biocloud1k_', 1:1000, '.r', stacked = TRUE)
bbs_footprint_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'footprint_', 1:100, '.r')
bbs_geoage_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'geoage_', 1:100, '.r')
bbs_soil_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'soil_', 1:100, '.r')
bbs_night_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'nightlight_', 1:100, '.r')
bbs_dhi_stats <- get_stats('/mnt/research/nasabio/data/bbs/geostats', 'dhi_', 1:100, '.r')

# Find and replace any NA entries in the list of data frames with a data frame of the same dimensions but filled with NAs

bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv', stringsAsFactors = FALSE)

bbs_elevation_stats <- replace_na_df(bbs_elevation_stats)
bbs_slope_stats <- replace_na_df(bbs_slope_stats)
bbs_tpi_stats <- replace_na_df(bbs_tpi_stats)
bbs_sin_aspect_stats <- replace_na_df(bbs_sin_aspect_stats)
bbs_cos_aspect_stats <- replace_na_df(bbs_cos_aspect_stats)
bbs_bioclim5k_stats <- replace_na_df(bbs_bioclim5k_stats)
bbs_bioclim1k_stats <- replace_na_df(bbs_bioclim1k_stats)
bbs_biocloud5k_stats <- replace_na_df(bbs_biocloud5k_stats)
bbs_biocloud1k_stats <- replace_na_df(bbs_biocloud1k_stats)
bbs_footprint_stats <- replace_na_df(bbs_footprint_stats)
bbs_geoage_stats <- replace_na_df(bbs_geoage_stats)
bbs_soil_stats <- replace_na_df(bbs_soil_stats)
bbs_night_stats <- replace_na_df(bbs_night_stats)
bbs_dhi_stats <- replace_na_df(bbs_dhi_stats)

# Replace variable names
bbs_elevation_stats <- replace_varname(bbs_elevation_stats, 'elevation')
bbs_slope_stats <- replace_varname(bbs_slope_stats, 'slope')
bbs_tpi_stats <- replace_varname(bbs_tpi_stats, 'TPI')
bbs_sin_aspect_stats <- replace_varname(bbs_sin_aspect_stats, 'sin_aspect')
bbs_cos_aspect_stats <- replace_varname(bbs_cos_aspect_stats, 'cos_aspect')
bbs_bioclim5k_stats <- replace_varname(bbs_bioclim5k_stats, paste0('bio', 1:19, '_5k'))
bbs_bioclim1k_stats <- replace_varname(bbs_bioclim1k_stats, paste0('bio', 1:19,'_1k'))
bbs_biocloud5k_stats <- replace_varname(bbs_biocloud5k_stats, rep(paste0('biocloud', 1:8, '_5k'), each = 8))
bbs_biocloud1k_stats <- replace_varname(bbs_biocloud1k_stats, rep(paste0('biocloud', 1:8, '_1k'), each = 8))
bbs_footprint_stats <- replace_varname(bbs_footprint_stats, 'human_footprint')
bbs_geoage_stats <- replace_varname(bbs_geoage_stats, 'geological_age')
bbs_soil_stats <- replace_varname(bbs_soil_stats, 'soil_type')
bbs_night_stats <- replace_varname(bbs_night_stats, 'nightlight')
bbs_dhi_stats <- replace_varname(bbs_dhi_stats, rep(c('dhi_fpar', 'dhi_gpp', 'dhi_lai8', 'dhi_ndvi'), each = 8))

bbs_all_stats <- list()

for (i in 1:length(bbs_elevation_stats)) {
	bbs_all_stats[[i]] <- full_join(cbind(bbsll[i,], rbind(bbs_elevation_stats[[i]],
														   bbs_slope_stats[[i]],
														   bbs_tpi_stats[[i]],
														   bbs_sin_aspect_stats[[i]],
														   bbs_cos_aspect_stats[[i]],
														   bbs_bioclim5k_stats[[i]],
														   bbs_bioclim1k_stats[[i]],
														   bbs_biocloud5k_stats[[i]],
														   bbs_biocloud1k_stats[[i]],
														   bbs_footprint_stats[[i]],
														   bbs_night_stats[[i]],
														   bbs_dhi_stats[[i]])),
									cbind(bbsll[i,], rbind(bbs_geoage_stats[[i]],
														   bbs_soil_stats[[i]])))
}

bbs_all_stats <- do.call('rbind', bbs_all_stats)
bbs_all_stats <- rename(bbs_all_stats, lon_aea = lon.1, lat_aea = lat.1, richness_geodiv = richness, diversity_geodiv = diversity)

write.csv(bbs_all_stats, file = '/mnt/research/nasabio/data/bbs/bbs_geodiversity_stats.csv', row.names = FALSE)

##############################
# FIA

fia_elev_stats <- get_all_elev_stats('/mnt/research/nasabio/data/fia/elevstats/big30m', 'stats_', 1:5000, '.r')
fia_elevation_stats <- fia_elev_stats$elev
fia_slope_stats <- fia_elev_stats$slope
fia_tpi_stats <- fia_elev_stats$TPI
fia_sin_aspect_stats <- fia_elev_stats$sin_aspect
fia_cos_aspect_stats <- fia_elev_stats$cos_aspect

fia_bioclim5k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'bioclim5k_', 1:1000, '.r')
fia_bioclim1k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'bioclim1k_', 1:5000, '.r')
fia_biocloud5k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'biocloud5k_', 1:1000, '.r', stacked = TRUE)
fia_biocloud1k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'biocloud1k_', 1:5000, '.r', stacked = TRUE)
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

fia_elevation_stats <- replace_na_df(fia_elevation_stats)
fia_slope_stats <- replace_na_df(fia_slope_stats)
fia_tpi_stats <- replace_na_df(fia_tpi_stats)
fia_sin_aspect_stats <- replace_na_df(fia_sin_aspect_stats)
fia_cos_aspect_stats <- replace_na_df(fia_cos_aspect_stats)
fia_bioclim5k_stats <- replace_na_df(fia_bioclim5k_stats)
fia_bioclim1k_stats <- replace_na_df(fia_bioclim1k_stats)
fia_biocloud5k_stats <- replace_na_df(fia_biocloud5k_stats)
fia_biocloud1k_stats <- replace_na_df(fia_biocloud1k_stats)
fia_footprint_stats <- replace_na_df(fia_footprint_stats)
fia_geoage_stats <- replace_na_df(fia_geoage_stats)
fia_soil_stats <- replace_na_df(fia_soil_stats)
fia_night_stats <- replace_na_df(fia_night_stats)
fia_dhi_stats <- replace_na_df(fia_dhi_stats)

# Replace variable names
fia_elevation_stats <- replace_varname(fia_elevation_stats, 'elevation')
fia_slope_stats <- replace_varname(fia_slope_stats, 'slope')
fia_tpi_stats <- replace_varname(fia_tpi_stats, 'TPI')
fia_sin_aspect_stats <- replace_varname(fia_sin_aspect_stats, 'sin_aspect')
fia_cos_aspect_stats <- replace_varname(fia_cos_aspect_stats, 'cos_aspect')
fia_bioclim5k_stats <- replace_varname(fia_bioclim5k_stats, paste0('bio', 1:19, '_5k'))
fia_bioclim1k_stats <- replace_varname(fia_bioclim1k_stats, paste0('bio', 1:19,'_1k'))
fia_biocloud5k_stats <- replace_varname(fia_biocloud5k_stats, rep(paste0('biocloud', 1:8, '_5k'), each = 8))
fia_biocloud1k_stats <- replace_varname(fia_biocloud1k_stats, rep(paste0('biocloud', 1:8, '_1k'), each = 8))
fia_footprint_stats <- replace_varname(fia_footprint_stats, 'human_footprint')
fia_geoage_stats <- replace_varname(fia_geoage_stats, 'geological_age')
fia_soil_stats <- replace_varname(fia_soil_stats, 'soil_type')
fia_night_stats <- replace_varname(fia_night_stats, 'nightlight')
fia_dhi_stats <- replace_varname(fia_dhi_stats, rep(c('dhi_fpar', 'dhi_gpp', 'dhi_lai8', 'dhi_ndvi'), each = 8))			

fia_all_stats <- list()

for (i in 1:length(fia_elevation_stats)) {
	fia_all_stats[[i]] <- full_join(cbind(fiall[i,], rbind(fia_elevation_stats[[i]],
														   fia_slope_stats[[i]],
														   fia_tpi_stats[[i]],
														   fia_sin_aspect_stats[[i]],
														   fia_cos_aspect_stats[[i]],
														   fia_bioclim5k_stats[[i]],
														   fia_bioclim1k_stats[[i]],
														   fia_biocloud5k_stats[[i]],
														   fia_biocloud1k_stats[[i]],
														   fia_footprint_stats[[i]],
														   fia_night_stats[[i]],
														   fia_dhi_stats[[i]])),
									cbind(fiall[i,], rbind(fia_geoage_stats[[i]],
														   fia_soil_stats[[i]])))
}


fia_all_stats <- do.call('rbind', fia_all_stats)
fia_all_stats <- rename(fia_all_stats, richness_geodiv = richness, diversity_geodiv = diversity)

write.csv(fia_all_stats, file = '/mnt/research/nasabio/data/fia/fia_geodiversity_stats.csv', row.names = FALSE)
