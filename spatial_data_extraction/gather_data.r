# Put together all diversity values for Andy's analysis in wide format
# QDR NASAbioxgeo 14 Feb 2018

library(dplyr)
library(reshape2)

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

fiageorad <- read.csv(file.path(fpdata, 'fia/fia_geodiversity_reduced.csv'), stringsAsFactors = FALSE)

fiageorad_wide <- fiageorad %>%
	filter(variable %in% c('bio1_5k', 'bio12_5k','elevation_5k', 'geological_age_5k')) %>%
	mutate(variable = gsub('_5k', '', variable, fixed = TRUE)) %>%
	select(PLT_CN, radius, mean, sd, variable, richness_geodiv, diversity_geodiv, mode) 

# Create wide format.
fiageorad_wide_sd <- fiageorad_wide %>% filter(!variable %in% 'geological_age') %>% dcast(PLT_CN ~ variable + radius, value.var = 'sd') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'sd', sep = '_')))
fiageorad_wide_mean <- fiageorad_wide %>% filter(!variable %in% 'geological_age') %>% dcast(PLT_CN ~ variable + radius, value.var = 'mean') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mean', sep = '_')))
fiageorad_wide_div <- fiageorad_wide %>% filter(variable %in% 'geological_age') %>% dcast(PLT_CN ~ variable + radius, value.var = 'diversity_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'diversity', sep = '_')))
fiageorad_wide_mode <- fiageorad_wide %>% filter(variable %in% 'geological_age') %>% dcast(PLT_CN ~ variable + radius, value.var = 'mode') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mode', sep = '_')))

# Biodiversity

# Point values
fiaalphapt <- read.csv(file.path(fpdata, 'fia/fiausa_alphadiv.csv'), stringsAsFactors = FALSE) %>% select(PLT_CN, shannon) %>% mutate(shannon = exp(shannon)) %>% rename(alpha_point = shannon)

# Radius values
fiaalpharadius <- read.csv(file.path(fpdata, 'fia/fiausa_alpha.csv'), stringsAsFactors = FALSE) %>% 
	select(PLT_CN, radius, shannon) %>%
	mutate(shannon = exp(shannon)) %>%
	rename(alpha = shannon) %>%
	dcast(PLT_CN ~ radius, value.var = 'alpha') %>%
	setNames(c(names(.)[1], paste('alpha', names(.)[-1], sep = '_')))

fiabetaradius <- read.csv(file.path(fpdata, 'fia/fiausa_betatd.csv'), stringsAsFactors = FALSE) %>% 
	select(PLT_CN, radius, beta_td_pairwise) %>% 
	rename(beta = beta_td_pairwise) %>%
	dcast(PLT_CN ~ radius, value.var = 'beta') %>%
	setNames(c(names(.)[1], paste('beta', names(.)[-1], sep = '_')))
	
fiagammaradius <- read.csv(file.path(fpdata, 'fia/fiausa_gamma.csv'), stringsAsFactors = FALSE) %>% 
	select(PLT_CN, radius, shannon) %>%
	mutate(shannon = exp(shannon)) %>%
	rename(gamma = shannon) %>%
	dcast(PLT_CN ~ radius, value.var = 'gamma') %>%
	setNames(c(names(.)[1], paste('gamma', names(.)[-1], sep = '_')))
	
# Join all

fiadat <- fiadat %>%
	left_join(fiageorad_wide_sd) %>%
	left_join(fiageorad_wide_mean) %>%
	left_join(fiageorad_wide_div) %>%
	left_join(fiageorad_wide_mode) %>%
	left_join(fiaalphapt) %>%
	left_join(fiaalpharadius) %>%
	left_join(fiabetaradius) %>%
	left_join(fiagammaradius)

names(fiadat) <- gsub('bio1_', 'MAT_', names(fiadat))	
names(fiadat) <- gsub('bio12_', 'MAP_', names(fiadat))
names(fiadat) <- gsub('geological_age', 'geoage', names(fiadat))		

write.csv(fiadat, file = file.path(fpdata, 'fia/fia_biogeo_joined.csv'), row.names = FALSE)

###########################################################

# BBS

# Geodiversity

# Point values

bbsgeopt <- read.csv(file.path(fpdata, 'bbs/bbs_geo_by_point.csv'), stringsAsFactors = FALSE)
bbshuc <- read.csv(file.path(fpdata, 'bbs/bbs_huc4.csv'), stringsAsFactors = FALSE)
bbsll <- read.csv(file.path(fpdata, 'bbs/bbs_correct_route_centroids.csv'), stringsAsFactors = FALSE)

bbsdat <- bbsgeopt %>%
	left_join(bbsll) %>%
	select(rteNo, lat, lon, elevation_30m, bio1_1k, bio12_1k, geological_age_1k) %>%
	rename(elevation_point = elevation_30m, MAT_point = bio1_1k, MAP_point = bio12_1k, geoage_point = geological_age_1k) %>%
	left_join(bbshuc)

# Radius values

bbsgeorad <- read.csv(file.path(fpdata, 'bbs/bbs_geodiversity_reduced.csv'), stringsAsFactors = FALSE)

bbsgeorad_wide <- bbsgeorad %>%
	filter(variable %in% c('bio1_5k', 'bio12_5k','elevation_5k', 'geological_age_5k')) %>%
	mutate(variable = gsub('_5k', '', variable, fixed = TRUE)) %>%
	select(rteNo, radius, mean, sd, variable, richness_geodiv, diversity_geodiv, mode) 

# Create wide format.
bbsgeorad_wide_sd <- bbsgeorad_wide %>% filter(!variable %in% 'geological_age') %>% dcast(rteNo ~ variable + radius, value.var = 'sd') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'sd', sep = '_')))
bbsgeorad_wide_mean <- bbsgeorad_wide %>% filter(!variable %in% 'geological_age') %>% dcast(rteNo ~ variable + radius, value.var = 'mean') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mean', sep = '_')))
bbsgeorad_wide_div <- bbsgeorad_wide %>% filter(variable %in% 'geological_age') %>% dcast(rteNo ~ variable + radius, value.var = 'diversity_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'diversity', sep = '_')))
bbsgeorad_wide_mode <- bbsgeorad_wide %>% filter(variable %in% 'geological_age') %>% dcast(rteNo ~ variable + radius, value.var = 'mode') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mode', sep = '_')))

# Biodiversity

# Point values
bbsalphapt <- read.csv(file.path(fpdata, 'bbs/bbs_alphadiv_1year.csv'), stringsAsFactors = FALSE) %>% select(rteNo, richness) %>% rename(alpha_point = richness)

# Radius values
bbsalpharadius <- read.csv(file.path(fpdata, 'bbs/bbs_alpha_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(rteNo, radius, richness) %>% 
	rename(alpha = richness) %>%
	dcast(rteNo ~ radius, value.var = 'alpha') %>%
	setNames(c(names(.)[1], paste('alpha', names(.)[-1], sep = '_')))

bbsbetaradius <- read.csv(file.path(fpdata, 'bbs/bbs_betatdpdfd_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(rteNo, radius, beta_td_pairwise_pa) %>% 
	rename(beta = beta_td_pairwise_pa) %>%
	dcast(rteNo ~ radius, value.var = 'beta') %>%
	setNames(c(names(.)[1], paste('beta', names(.)[-1], sep = '_')))
	
bbsgammaradius <- read.csv(file.path(fpdata, 'bbs/bbs_gamma_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(rteNo, radius, richness) %>% 
	rename(gamma = richness) %>%
	dcast(rteNo ~ radius, value.var = 'gamma') %>%
	setNames(c(names(.)[1], paste('gamma', names(.)[-1], sep = '_')))
	
# Join all

bbsdat <- bbsdat %>%
	left_join(bbsgeorad_wide_sd) %>%
	left_join(bbsgeorad_wide_mean) %>%
	left_join(bbsgeorad_wide_div) %>%
	left_join(bbsgeorad_wide_mode) %>%
	left_join(bbsalphapt) %>%
	left_join(bbsalpharadius) %>%
	left_join(bbsbetaradius) %>%
	left_join(bbsgammaradius)

names(bbsdat) <- gsub('bio1_', 'MAT_', names(bbsdat))	
names(bbsdat) <- gsub('bio12_', 'MAP_', names(bbsdat))
names(bbsdat) <- gsub('geological_age', 'geoage', names(bbsdat))		

write.csv(bbsdat, file = file.path(fpdata, 'bbs/bbs_biogeo_joined.csv'), row.names = FALSE)

###############################################################

# 21 Feb 2018
# All geodiversity variables in a single wide format data frame for BBS.

library(dplyr)
library(reshape2)
fpdata <- '/mnt/research/nasabio/data'

bbsgeopt <- read.csv(file.path(fpdata, 'bbs/bbs_geo_by_point.csv'), stringsAsFactors = FALSE)
bbshuc <- read.csv(file.path(fpdata, 'bbs/bbs_huc4.csv'), stringsAsFactors = FALSE)
bbsll <- read.csv(file.path(fpdata, 'bbs/bbs_correct_route_centroids.csv'), stringsAsFactors = FALSE)

bbsdat <- bbsgeopt %>%
	left_join(bbsll) %>%
	left_join(bbshuc) %>%
	select(-lon.1, -lat.1) %>%
	select(rteNo, lat, lon, HUC4, everything()) %>%
	setNames(c(names(.)[1:4], paste(names(.)[-(1:4)], 'point', sep = '_')))

# Radius values

bbsgeorad <- read.csv(file.path(fpdata, 'bbs/bbs_geodiversity.csv'), stringsAsFactors = FALSE)

bbsgeorad <- bbsgeorad %>%
	select(rteNo, radius, mean, sd, variable, richness_geodiv, diversity_geodiv, mode) 

discrete_vars <- grep('soil|geol', unique(bbsgeorad$variable), value = TRUE)
	
# Create wide format.
bbsgeorad_wide_sd <- bbsgeorad %>% filter(!variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'sd') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'sd', sep = '_')))
bbsgeorad_wide_mean <- bbsgeorad %>% filter(!variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'mean') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mean', sep = '_')))
bbsgeorad_wide_div <- bbsgeorad %>% filter(variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'diversity_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'diversity', sep = '_')))
bbsgeorad_wide_rich <- bbsgeorad %>% filter(variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'richness_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'richness', sep = '_')))
bbsgeorad_wide_mode <- bbsgeorad %>% filter(variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'mode') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mode', sep = '_')))

# Join all

bbsdat <- bbsdat %>%
	left_join(bbsgeorad_wide_sd) %>%
	left_join(bbsgeorad_wide_mean) %>%
	left_join(bbsgeorad_wide_div) %>%
	left_join(bbsgeorad_wide_rich) %>%
	left_join(bbsgeorad_wide_mode)
	
write.csv(bbsdat, file = file.path(fpdata, 'bbs/bbs_allgeo_wide.csv'), row.names = FALSE)

# All biodiversity in a single wide format data frame for BBS.
bbsalphapt <- read.csv(file.path(fpdata, 'bbs/bbs_alphadiv_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(rteNo, lon, lat, lon_aea, lat_aea, richness, MPD_pa_z, MNTD_pa_z, MPDfunc_pa_z, MNTDfunc_pa_z) %>%
	setNames(c(names(.)[1:5], paste('alpha', names(.)[-(1:5)], 'point', sep = '_')))

library(data.table)

# Radius values
bbsalpharadius <- read.csv(file.path(fpdata, 'bbs/bbs_alpha_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(-lon, -lat, -lon_aea, -lat_aea) %>% 
	setDT %>%
	dcast(rteNo ~ radius, value.var = c("richness", "MPD_pa_z", "MNTD_pa_z", "MPDfunc_pa_z", "MNTDfunc_pa_z")) %>%
	setNames(c(names(.)[1], paste('alpha', names(.)[-1], sep = '_')))

bbsbetaradius <- read.csv(file.path(fpdata, 'bbs/bbs_betatdpdfd_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(-lon, -lat, -lon_aea, -lat_aea) %>% 
	setDT %>%
	dcast(rteNo ~ radius, value.var = c("beta_td_pairwise_pa", "beta_td_sorensen_pa",  "beta_pd_pairwise_pa", "beta_pd_pairwise_pa_z", "beta_pd_nt_pa", "beta_pd_nt_pa_z", "beta_fd_pairwise_pa", "beta_fd_pairwise_pa_z", "beta_fd_nt_pa", "beta_fd_nt_pa_z"))
	
bbsgammaradius <- read.csv(file.path(fpdata, 'bbs/bbs_gamma_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(-lon, -lat, -lon_aea, -lat_aea) %>% 
	setDT %>%
	dcast(rteNo ~ radius, value.var = c("richness", "MPD_pa_z", "MNTD_pa_z", "MPDfunc_pa_z", "MNTDfunc_pa_z")) %>%
	setNames(c(names(.)[1], paste('gamma', names(.)[-1], sep = '_')))
	
bbsbiodat <- bbsalphapt %>%
	left_join(bbsalpharadius) %>%
	left_join(bbsbetaradius) %>%
	left_join(bbsgammaradius)
	
write.csv(bbsbiodat, file = file.path(fpdata, 'bbs/bbs_allbio_wide.csv'), row.names = FALSE)