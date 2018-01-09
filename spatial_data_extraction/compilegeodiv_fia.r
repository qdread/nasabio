# FIA compile geodiversity stats (entire USA)
# Submit job because it is too big to run on dev node, as usual

# Edited 09 Jan: compile elevation only for entire USA (ignore other variables for the moment)

library(dplyr)

# Define functions

get_stats <- function(fp, prefix, n, suffix, stacked = FALSE) {
	file_names <- paste0(prefix, n, suffix)
	all_stats <- list()
	pb <- txtProgressBar(0, length(file_names), style = 3)
	for (i in 1:length(file_names)) {
		setTxtProgressBar(pb, i)
		load(file.path(fp, file_names[i]))
		if (!stacked) all_stats[[length(all_stats) + 1]] <- stats_by_point
		if (stacked) all_stats[[length(all_stats) + 1]] <- lapply(stats_by_point, function(x) do.call('rbind', x))
	}
	close(pb)
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

fia_path <- '/mnt/research/nasabio/data/fia/elevstats_usa'

fia_elevation_stats <- get_stats(fia_path, 'elevation_', 1:10000, '.r', stacked = TRUE)
fia_slope_stats <- get_stats(fia_path, 'slope_', 1:10000, '.r', stacked = TRUE)
fia_tpi_stats <- get_stats(fia_path, 'tri_', 1:10000, '.r', stacked = TRUE)
fia_aspect_stats <- get_stats(fia_path, 'roughness_', 1:10000, '.r', stacked = TRUE)

# Add aspect here.

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
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa.csv', stringsAsFactors = FALSE)
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
fia_biocloud5k_stats <- replace_varname(fia_biocloud5k_stats, paste0('biocloud', 1:8, '_5k'))
fia_biocloud1k_stats <- replace_varname(fia_biocloud1k_stats, paste0('biocloud', 1:8, '_1k'))
fia_footprint_stats <- replace_varname(fia_footprint_stats, 'human_footprint')
fia_geoage_stats <- replace_varname(fia_geoage_stats, 'geological_age')
fia_soil_stats <- replace_varname(fia_soil_stats, 'soil_type')
fia_night_stats <- replace_varname(fia_night_stats, 'nightlight')
fia_dhi_stats <- replace_varname(fia_dhi_stats, c('dhi_fpar', 'dhi_gpp', 'dhi_lai8', 'dhi_ndvi'))			

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


# Edit 11 Dec. Use the much faster dplyr analog to do.call rbind
fia_all_stats <- bind_rows(fia_all_stats)

fia_all_stats <- rename(fia_all_stats, richness_geodiv = richness, diversity_geodiv = diversity)

write.csv(fia_all_stats, file = '/mnt/research/nasabio/data/fia/geodiv/fia_pnw_geodiversity_stats.csv', row.names = FALSE)

# Edit 11 Dec. Split this up into a few groups so that only needed variables can be loaded at one time.
elev_vars <- c('elevation','slope','TPI','sin_aspect','cos_aspect')
bio5k_vars <- paste0('bio', 1:19, '_5k')
bio1k_vars <- paste0('bio', 1:19, '_1k')
biocloud_vars <- c(paste0('biocloud', 1:8, '_5k'), paste0('biocloud', 1:8, '_1k'))
other_vars <- c('human_footprint', 'geological_age', 'soil_type', 'nightlight', 'dhi_fpar', 'dhi_gpp', 'dhi_lai8', 'dhi_ndvi')

write.csv(filter(fia_all_stats, variable %in% elev_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_pnw_elev_stats.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% bio5k_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_pnw_bio5k_stats.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% bio1k_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_pnw_bio1k_stats.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% biocloud_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_pnw_biocloud_stats.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% other_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_pnw_other_stats.csv', row.names = FALSE)

# Added 09 Jan: write elevation only.
fia_elevation_stats <- cbind(PLT_CN = rep(plotmetadata$PLT_CN, each = 11), bind_rows(fia_elevation_stats))
write.csv(fia_elevation_stats, file = '/mnt/research/nasabio/data/fia/geodiv/fia_usa_elev_only.csv', row.names = FALSE)