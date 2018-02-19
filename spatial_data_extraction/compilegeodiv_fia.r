# FIA compile geodiversity stats (entire USA)
# Submit job because it is too big to run on dev node, as usual

# Edited 19 Feb: use variable table to get all variable names and number of files for each taxon
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

vartable <- read.csv('/mnt/research/nasabio/data/geodiv_table_for_gdal_reorder.csv', stringsAsFactors = FALSE)
correct_var_names <- list('elevation_30m', 'elevation_30m_tri', 'elevation_30m_roughness',
						  'elevation_5k', 'elevation_5k_tri', 'elevation_5k_roughness',
						  'slope_30m', 'slope_5k',
						  'aspect_sin_30m', 'aspect_cos_30m', 'aspect_sin_5k', 'aspect_cos_5k',
						  paste0('bio', 1:19,'_1k'),
						  paste0('bio', 1:19,'_1k_tri'),
						  paste0('bio', 1:19,'_1k_roughness'),
						  paste0('bio', 1:19,'_5k'),
						  paste0('bio', 1:19,'_5k_tri'),
						  paste0('bio', 1:19,'_5k_roughness'),
						  paste0('biocloud', 1:8,'_1k'),
						  paste0('biocloud', 1:8,'_1k_tri'),
						  paste0('biocloud', 1:8,'_1k_roughness'),
						  paste0('biocloud', 1:8,'_5k'),
						  paste0('biocloud', 1:8,'_5k_tri'),
						  paste0('biocloud', 1:8,'_5k_roughness'),
						  c('dhi_fpar_1k', 'dhi_gpp_1k', 'dhi_lai8_1k', 'dhi_ndvi_1k'),
						  c('dhi_fpar_1k_tri', 'dhi_gpp_1k_tri', 'dhi_lai8_1k_tri', 'dhi_ndvi_1k_tri'),
						  c('dhi_fpar_1k_roughness', 'dhi_gpp_1k_roughness', 'dhi_lai8_1k_roughness', 'dhi_ndvi_1k_roughness'),
						  c('dhi_fpar_5k', 'dhi_gpp_5k', 'dhi_lai8_5k', 'dhi_ndvi_5k'),
						  c('dhi_fpar_5k_tri', 'dhi_gpp_5k_tri', 'dhi_lai8_5k_tri', 'dhi_ndvi_5k_tri'),
						  c('dhi_fpar_5k_roughness', 'dhi_gpp_5k_roughness', 'dhi_lai8_5k_roughness', 'dhi_ndvi_5k_roughness'),
						  'human_footprint_1k', 'human_footprint_1k_tri', 'human_footprint_1k_roughness',
						  'human_footprint_5k', 'human_footprint_5k_tri', 'human_footprint_5k_roughness',
						  'nightlight_500m', 'nighlight_500m_tri', 'nightlight_500m_roughness',
						  'nightlight_5k', 'nightlight_5k_tri', 'nightlight_5k_roughness',
						  'geological_age_1k', 'geological_age_5k', 'soil_type_5k')

# READ BBS DATA #					
						  
bbs_path <- '/mnt/research/nasabio/data/bbs/allgeodiv_v2'
bbs_stats <- list()

for (i in 1:nrow(vartable)) {
	bbs_stats[[i]] <- get_stats(bbs_path, paste0(vartable$variable.id[i], '_'), 1:vartable$N.slices.bbs, '.r', stacked = TRUE)
	bbs_stats[[i]] <- replace_na_df(bbs_stats[[i]])
	bbs_stats[[i]] <- replace_varname(bbs_stats[[i]], correct_var_names[[i]])
}

bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv', stringsAsFactors = FALSE)

for (i in 1:length(bbs_stats[[1]])) {
	stats_i <- lapply(bbs_stats, '[[', i)
	bbs_all_stats[[i]] <- full_join(cbind(plotmetadata[i,], bind_rows(stats_i[1:42])),
									cbind(plotmetadata[i,], bind_rows(stats_i[43:45])))
}

bbs_all_stats <- bind_rows(bbs_all_stats)
bbs_all_stats <- rename(bbs_all_stats, richness_geodiv = richness, diversity_geodiv = diversity)

write.csv(bbs_all_stats, file = '/mnt/research/nasabio/data/bbs/bbs_geodiversity.csv', row.names = FALSE)

# READ FIA DATA #

fia_path <- '/mnt/research/nasabio/data/fia/allgeodiv_v2'
fia_stats <- list()


for (i in 1:nrow(vartable)) {
	fia_stats[[i]] <- get_stats(fia_path, paste0(vartable$variable.id[i], '_'), 1:vartable$N.slices.fiaall[i], '.r', stacked = TRUE)
	fia_stats[[i]] <- replace_na_df(fia_stats[[i]])
	fia_stats[[i]] <- replace_varname(fia_stats[[i]], correct_var_names[[i]])
}

plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa.csv', stringsAsFactors = FALSE)

fia_all_stats <- list()

for (i in 1:length(fia_stats[[1]])) {
	stats_i <- lapply(fia_stats, '[[', i)
	fia_all_stats[[i]] <- full_join(cbind(plotmetadata[i,], bind_rows(stats_i[1:42])),
									cbind(plotmetadata[i,], bind_rows(stats_i[43:45])))
}

fia_all_stats <- bind_rows(fia_all_stats)
fia_all_stats <- rename(fia_all_stats, richness_geodiv = richness, diversity_geodiv = diversity)

write.csv(fia_all_stats, file = '/mnt/research/nasabio/data/fia/geodiv/fia_usa_geodiversity.csv', row.names = FALSE)

# Also split it up so that the files are reasonably sized.
elev_vars <- unlist(correct_var_names[1:12])
bio1k_vars<- unlist(correct_var_names[13:15])
bio5k_vars <- unlist(correct_var_names[16:18])
biocloud_vars <- unlist(correct_var_names[19:24])
other_vars <- unlist(correct_var_names[25:45])

write.csv(filter(fia_all_stats, variable %in% elev_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_usa_elev.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% bio5k_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_usa_bio5k.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% bio1k_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_usa_bio1k.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% biocloud_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_usa_biocloud.csv', row.names = FALSE)
write.csv(filter(fia_all_stats, variable %in% other_vars), file = '/mnt/research/nasabio/data/fia/geodiv/fia_usa_other.csv', row.names = FALSE)

