# Load and combine 1-km geodiversity values
# QDR / NASAbioXgeo / 1 July 2018

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

correct_var_names <- list('aspect_cos_30m', 'aspect_sin_30m', 
                          paste0('bio', 1:19, '_1k'),
                          paste0('bio', 1:19,'_1k_roughness'),
                          paste0('bio', 1:19,'_1k_tri'),
                          paste0('biocloud', 1:8,'_1k'),
                          paste0('biocloud', 1:8,'_1k_roughness'),
                          paste0('biocloud', 1:8,'_1k_rtri'),
                          c('dhi_fpar_1k', 'dhi_gpp_1k', 'dhi_lai8_1k', 'dhi_ndvi_1k'),
                          c('dhi_fpar_1k_roughness', 'dhi_gpp_1k_roughness', 'dhi_lai8_1k_roughness', 'dhi_ndvi_1k_roughness'),
                          c('dhi_fpar_1k_tri', 'dhi_gpp_1k_tri', 'dhi_lai8_1k_tri', 'dhi_ndvi_1k_tri'),
                          'elevation_30m', 
                          'geological_age_1k',
                          'human_footprint_1k', 'human_footprint_1k_roughness', 'human_footprint_1k_tri',
                          'nightlight_500m', 'nightlight_500m_roughness', 'nightlight_500m_tri',
                          'elevation_30m_roughness', 'slope_30m', 'elevation_30m_tri')

var_names <- c("aspect_cos", "aspect_sin", "bioclim1k", "bioclim1k_roughness",
               "bioclim1k_tri", "biocloud1k", "biocloud1k_roughness", "biocloud1k_tri",
               "dhi", "dhi_roughness", "dhi_tri", "elevation", "gea", "hf",
               "hf_roughness", "hf_tri", "night", "night_roughness", "night_tri",
               "roughness", "slope", "tri")

# READ BBS DATA #					

bbs_path <- '/mnt/research/nasabio/data/bbs/allgeodiv_1km'
bbs_stats <- list()

for (i in 1:length(var_names)) {
  bbs_stats[[i]] <- get_stats(bbs_path, paste0(var_names[i], '_'), 1, '.r', stacked = TRUE)
  bbs_stats[[i]] <- replace_na_df(bbs_stats[[i]])
  bbs_stats[[i]] <- replace_varname(bbs_stats[[i]], correct_var_names[[i]])
}

bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv', stringsAsFactors = FALSE)

bbs_all_stats <- list()

for (i in 1:length(bbs_stats[[1]])) {
  stats_i <- lapply(bbs_stats, '[[', i)
  bbs_all_stats[[i]] <- full_join(cbind(bbsll[i,], bind_rows(stats_i[c(1:12, 14:22)])),
                                  cbind(bbsll[i,], bind_rows(stats_i[13])))
}

bbs_all_stats <- bind_rows(bbs_all_stats)
bbs_all_stats <- bbs_all_stats %>%
  rename(richness_geodiv = richness, diversity_geodiv = diversity) %>%
  select(-layer) %>%
  select(rteNo, lon, lat, lon_aea, lat_aea, variable, radius, everything())

write.csv(bbs_all_stats, file = '/mnt/research/nasabio/data/bbs/geodiversity_CSVs/bbs_geodiversity_1kmradius.csv', row.names = FALSE)

# READ FIA DATA #					

fia_path <- '/mnt/research/nasabio/data/fia/allgeodiv_1km'
fia_stats <- list()

for (i in 1:length(var_names)) {
  fia_stats[[i]] <- get_stats(fia_path, paste0(var_names[i], '_'), 1:50, '.r', stacked = TRUE)
  fia_stats[[i]] <- replace_na_df(fia_stats[[i]])
  fia_stats[[i]] <- replace_varname(fia_stats[[i]], correct_var_names[[i]])
}

plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa.csv', stringsAsFactors = FALSE)

fia_all_stats <- list()

for (i in 1:length(fia_stats[[1]])) {
  stats_i <- lapply(fia_stats, '[[', i)
  fia_all_stats[[i]] <- full_join(cbind(PLT_CN = plotmetadata[i,], bind_rows(stats_i[c(1:12,14:22)])),
                                  cbind(PLT_CN = plotmetadata[i,], bind_rows(stats_i[13])))
}

fia_all_stats <- bind_rows(fia_all_stats)
fia_all_stats <- fia_all_stats %>%
  rename(richness_geodiv = richness, diversity_geodiv = diversity) %>%
  select(-layer) %>%
  select(PLT_CN, variable, radius, everything())

write.csv(fia_all_stats, file = '/mnt/research/nasabio/data/fia/geodiversity_CSVs/fia_geodiversity_1kmradius.csv', row.names = FALSE)