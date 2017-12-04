# Compile summary statistics: elevation and bioclim.
# BBS and FIA
# 04 Sep 2017 QDR

# Edited 01 Dec 2017: updated for new workflow (just for FIA elevation, other datasets will be updated as needed)
# Edited 04 Oct 2017: add new elevation stats that include elevation, slope, sin aspect, cos aspect, and TPI
# Edited 25 Sep 2017: add biocloud
# Edited 13 Sep 2017: add geostats (footprint, soil type diversity, and geological age diversity)

library(dplyr)

# Define functions

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

##############################
# BBS

bbs_path <- '/mnt/research/nasabio/data/bbs/allgeodiv'

bbs_elevation_stats <- get_stats(bbs_path, 'elevation_', 1:2000, '.r')
bbs_slope_stats <- get_stats(bbs_path, 'slope_', 1:2000, '.r')
bbs_tpi_stats <- get_stats(bbs_path, 'tpi_', 1:2000, '.r')
bbs_aspect_stats <- get_stats(bbs_path, 'aspect_', 1:2000, '.r')

bbs_bioclim5k_stats <- get_stats(bbs_path, 'bioclim5k_', 1:500, '.r')
bbs_bioclim1k_stats <- get_stats(bbs_path, 'bioclim1k_', 1:1000, '.r')
bbs_biocloud5k_stats <- get_stats(bbs_path, 'biocloud5k_', 1:500, '.r')
bbs_biocloud1k_stats <- get_stats(bbs_path, 'biocloud1k_', 1:1000, '.r')
bbs_footprint_stats <- get_stats(bbs_path, 'hf_', 1:100, '.r')
bbs_geoage_stats <- get_stats(bbs_path, 'gea_', 1:100, '.r')
bbs_soil_stats <- get_stats(bbs_path, 'soil_', 1:100, '.r')
bbs_night_stats <- get_stats(bbs_path, 'night_', 1:100, '.r')
bbs_dhi_stats <- get_stats(bbs_path, 'dhi_', 1:100, '.r')

# Find and replace any NA entries in the list of data frames with a data frame of the same dimensions but filled with NAs

bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv', stringsAsFactors = FALSE)

bbs_elevation_stats <- replace_na_df(bbs_elevation_stats)
bbs_slope_stats <- replace_na_df(bbs_slope_stats)
bbs_tpi_stats <- replace_na_df(bbs_tpi_stats)
bbs_aspect_stats <- replace_na_df(bbs_aspect_stats)
bbs_bioclim5k_stats <- replace_na_df(bbs_bioclim5k_stats)
bbs_bioclim1k_stats <- replace_na_df(bbs_bioclim1k_stats)
bbs_biocloud5k_stats <- replace_na_df(bbs_biocloud5k_stats)
bbs_biocloud1k_stats <- replace_na_df(bbs_biocloud1k_stats)
bbs_footprint_stats <- replace_na_df(bbs_footprint_stats)
bbs_geoage_stats <- replace_na_df(bbs_geoage_stats)
bbs_soil_stats <- replace_na_df(bbs_soil_stats)
bbs_night_stats <- replace_na_df(bbs_night_stats)
bbs_dhi_stats <- replace_na_df(bbs_dhi_stats)

bbs_sin_aspect_stats <- lapply(bbs_aspect_stats, function(x) with(x, data.frame(radius=radius, variable=variable, mean=mean_sin, sd=sd_sin, min=min_sin, max=max_sin, n=n_sin)))
bbs_cos_aspect_stats <- lapply(bbs_aspect_stats, function(x) with(x, data.frame(radius=radius, variable=variable, mean=mean_cos, sd=sd_cos, min=min_cos, max=max_cos, n=n_cos)))


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

fia_path <- '/mnt/research/nasabio/data/fia/allgeodiv'

fia_elevation_stats <- get_stats(fia_path, 'elevation_', 1:22531, '.r')
fia_slope_stats <- get_stats(fia_path, 'slope_', 1:22531, '.r')
fia_tpi_stats <- get_stats(fia_path, 'tpi_', 1:22531, '.r')
fia_aspect_stats <- get_stats(fia_path, 'aspect_', 1:22531, '.r')

fia_bioclim5k_stats <- get_stats(fia_path, 'bioclim5k_', 1:1000, '.r')
fia_bioclim1k_stats <- get_stats(fia_path, 'bioclim1k_', 1:5000, '.r')
fia_biocloud5k_stats <- get_stats(fia_path, 'biocloud5k_', 1:1000, '.r')
fia_biocloud1k_stats <- get_stats(fia_path, 'biocloud1k_', 1:5000, '.r')
fia_footprint_stats <- get_stats(fia_path, 'hf_', 1:250, '.r')
fia_geoage_stats <- get_stats(fia_path, 'gea_', 1:250, '.r')
fia_soil_stats <- get_stats(fia_path, 'soil_', 1:250, '.r')
fia_night_stats <- get_stats(fia_path, 'night_', 1:250, '.r')
fia_dhi_stats <- get_stats(fia_path, 'dhi_', 1:250, '.r')

# FIA plot IDs (no coords)
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords.csv', stringsAsFactors = FALSE)
#source('/mnt/research/nasabio/code/loadfia.r')

fia_elevation_stats <- replace_na_df(fia_elevation_stats)
fia_slope_stats <- replace_na_df(fia_slope_stats)
fia_tpi_stats <- replace_na_df(fia_tpi_stats)
fia_aspect_stats <- replace_na_df(fia_aspect_stats)
fia_bioclim5k_stats <- replace_na_df(fia_bioclim5k_stats)
fia_bioclim1k_stats <- replace_na_df(fia_bioclim1k_stats)
fia_biocloud5k_stats <- replace_na_df(fia_biocloud5k_stats)
fia_biocloud1k_stats <- replace_na_df(fia_biocloud1k_stats)
fia_footprint_stats <- replace_na_df(fia_footprint_stats)
fia_geoage_stats <- replace_na_df(fia_geoage_stats)
fia_soil_stats <- replace_na_df(fia_soil_stats)
fia_night_stats <- replace_na_df(fia_night_stats)
fia_dhi_stats <- replace_na_df(fia_dhi_stats)

fia_sin_aspect_stats <- lapply(fia_aspect_stats, function(x) with(x, data.frame(radius=radius, variable=variable, mean=mean_sin, sd=sd_sin, min=min_sin, max=max_sin, n=n_sin)))
fia_cos_aspect_stats <- lapply(fia_aspect_stats, function(x) with(x, data.frame(radius=radius, variable=variable, mean=mean_cos, sd=sd_cos, min=min_cos, max=max_cos, n=n_cos)))

# Replace variable names
fia_elevation_stats <- replace_varname(fia_elevation_stats, 'elevation')
fia_slope_stats <- replace_varname(fia_slope_stats, 'slope')
fia_tpi_stats <- replace_varname(fia_tpi_stats, 'TPI')
fia_sin_aspect_stats <- replace_varname(fia_sin_aspect_stats, 'sin_aspect')
fia_cos_aspect_stats <- replace_varname(fia_cos_aspect_stats, 'cos_aspect')
fia_bioclim5k_stats <- replace_varname(fia_bioclim5k_stats, paste0('bio', 1:19, '_5k'))
fia_bioclim1k_stats <- replace_varname(fia_bioclim1k_stats, paste0('bio', 1:19,'_1k'))
fia_biocloud5k_stats <- replace_varname(fia_biocloud5k_stats, rep(paste0('biocloud', 1:8, '_5k'), each = 13))
fia_biocloud1k_stats <- replace_varname(fia_biocloud1k_stats, rep(paste0('biocloud', 1:8, '_1k'), each = 13))
fia_footprint_stats <- replace_varname(fia_footprint_stats, 'human_footprint')
fia_geoage_stats <- replace_varname(fia_geoage_stats, 'geological_age')
fia_soil_stats <- replace_varname(fia_soil_stats, 'soil_type')
fia_night_stats <- replace_varname(fia_night_stats, 'nightlight')
fia_dhi_stats <- replace_varname(fia_dhi_stats, rep(c('dhi_fpar', 'dhi_gpp', 'dhi_lai8', 'dhi_ndvi'), each = 13))			

fia_all_stats <- list()

for (i in 1:length(fia_elevation_stats)) {
	fia_all_stats[[i]] <- full_join(cbind(plotmetadata[i,], rbind(fia_elevation_stats[[i]],
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
									cbind(plotmetadata[i,], rbind(fia_geoage_stats[[i]],
														   fia_soil_stats[[i]])))
}


fia_all_stats <- do.call('rbind', fia_all_stats)
fia_all_stats <- rename(fia_all_stats, richness_geodiv = richness, diversity_geodiv = diversity)

write.csv(fia_all_stats, file = '/mnt/research/nasabio/data/fia/fia_geodiversity_stats.csv', row.names = FALSE)

# Temporary: just DEM (01 Dec)
fia_elevation_stats <- get_stats(fia_path, 'elevation_', 1:22531, '.r')
fia_elevation_stats <- replace_na_df(fia_elevation_stats)
fia_elevation_stats <- replace_varname(fia_elevation_stats, 'elevation')
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords.csv', stringsAsFactors = FALSE)

for (i in 1:length(fia_elevation_stats)) {
	fia_elevation_stats[[i]] <- cbind(plotmetadata[rep(i, nrow(fia_elevation_stats[[i]])),], fia_elevation_stats[[i]])
}
fia_elevation_stats <- do.call('rbind', fia_elevation_stats)
write.csv(fia_elevation_stats, file = '/mnt/research/nasabio/data/fia/fia_elev_stats_unfuzzed.csv', row.names = FALSE)
