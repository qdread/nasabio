
# 21 Feb 2018
# All geodiversity variables in a single wide format data frame for BBS and FIA.

# Edit 25 Apr 2018: Use midpoints, not centroids
# Edit 01 Jul 2018: add the 1 kilometer ones too

library(dplyr)
library(reshape2)
fpdata <- '/mnt/research/nasabio/data'

bbsgeopt <- read.csv(file.path(fpdata, 'bbs/bbs_geo_by_point.csv'), stringsAsFactors = FALSE)
bbseco <- read.csv(file.path(fpdata, 'bbs/bbs_ecoregions.csv'), stringsAsFactors = FALSE)
bbsll <- read.csv(file.path(fpdata, 'bbs/bbs_route_midpoints.csv'), stringsAsFactors = FALSE)

bbsdat <- bbsgeopt %>%
	left_join(bbsll) %>%
	left_join(bbseco) %>%
	select(rteNo, lat, lon, lat_aea, lon_aea, HUC4, BCR, TNC, everything()) %>%
	setNames(c(names(.)[1:8], paste(names(.)[-(1:8)], 'point', sep = '_')))

# Radius values

bbsgeorad <- read.csv(file.path(fpdata, 'bbs/bbs_geodiversity.csv'), stringsAsFactors = FALSE) # for all but 1 km

bbsgeorad <- bbsgeorad %>%
	select(rteNo, radius, mean, sd, variable, richness_geodiv, diversity_geodiv, mode) 

discrete_vars <- grep('soil|geol', unique(bbsgeorad$variable), value = TRUE)
tri_vars <- grep('_tri|_roughness', unique(bbsgeorad$variable), value = TRUE)
	
# Create wide format.
bbsgeorad_wide_sd <- bbsgeorad %>% filter(!variable %in% c(discrete_vars, tri_vars)) %>% dcast(rteNo ~ variable + radius, value.var = 'sd') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'sd', sep = '_')))
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
	left_join(bbsgeorad_wide_mode) %>%

# correct typo
bbsdat <- setNames(bbsdat, gsub('nighlight', 'nightlight', names(bbsdat)))
	
write.csv(bbsdat, file = file.path(fpdata, 'bbs/bbs_allgeo_wide.csv'), row.names = FALSE) # all but 1 km

# BBS geodiversity 1 km edition
bbsgeorad <- read.csv(file.path(fpdata, 'bbs/geodiversity_CSVs/bbs_geodiversity_1kmradius.csv'), stringsAsFactors = FALSE) # for 1 km

bbsgeorad <- bbsgeorad %>%
	select(rteNo, radius, mean, sd, variable, richness_geodiv, diversity_geodiv, mode) 

discrete_vars <- grep('soil|geol', unique(bbsgeorad$variable), value = TRUE)
tri_vars <- grep('_tri|_roughness', unique(bbsgeorad$variable), value = TRUE)
	
# Create wide format.
bbsgeorad_wide_sd <- bbsgeorad %>% filter(!variable %in% c(discrete_vars, tri_vars)) %>% dcast(rteNo ~ variable + radius, value.var = 'sd') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'sd', sep = '_')))
bbsgeorad_wide_mean <- bbsgeorad %>% filter(!variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'mean') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mean', sep = '_')))
bbsgeorad_wide_div <- bbsgeorad %>% filter(variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'diversity_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'diversity', sep = '_')))
bbsgeorad_wide_rich <- bbsgeorad %>% filter(variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'richness_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'richness', sep = '_')))
bbsgeorad_wide_mode <- bbsgeorad %>% filter(variable %in% discrete_vars) %>% dcast(rteNo ~ variable + radius, value.var = 'mode') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mode', sep = '_')))

# Join all

bbsdat <- bbsgeorad_wide_sd %>%
	left_join(bbsgeorad_wide_mean) %>%
	left_join(bbsgeorad_wide_div) %>%
	left_join(bbsgeorad_wide_rich) %>%
	left_join(bbsgeorad_wide_mode) %>%
  setNames(gsub('_rtri_', '_tri_', names(.))) # correct typo

write.csv(bbsdat, file = file.path(fpdata, 'bbs/geodiversity_CSVs/bbs_allgeo_wide_1kmradius.csv'), row.names = FALSE)


# All biodiversity in a single wide format data frame for BBS.
bbsalphapt <- read.csv(file.path(fpdata, 'bbs/bbs_alphadiv_1year.csv'), stringsAsFactors = FALSE) %>% 
	select(rteNo, richness, MPD_pa_z, MNTD_pa_z, MPDfunc_pa_z, MNTDfunc_pa_z) %>%
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

#############################################################
# This section added 24 Apr. 2018
# Compile within-route biodiversity for BBS
# Only includes radii 5, 10, and 20 km.

library(dplyr)
library(reshape2)
library(data.table)
fpdata <- '/mnt/research/nasabio/data'

bbsalphawithin <- read.csv(file.path(fpdata, 'bbs/biodiversity_CSVs/bbs_withinroute_alpha.csv'), stringsAsFactors = FALSE) %>%
	setDT %>%
	dcast(rteNo ~ radius, value.var = c("richness", "MPD_pa_z", "MNTD_pa_z", "MPDfunc_pa_z", "MNTDfunc_pa_z")) %>%
	setNames(c(names(.)[1], paste('alpha', names(.)[-1], sep = '_')))

bbsbetawithin <- read.csv(file.path(fpdata, 'bbs/biodiversity_CSVs/bbs_withinroute_beta.csv'), stringsAsFactors = FALSE) %>%
	setDT %>%
	dcast(rteNo ~ radius, value.var = c("beta_td_pairwise_pa", "beta_td_sorensen_pa",  "beta_pd_pairwise_pa", "beta_pd_pairwise_pa_z", "beta_pd_nt_pa", "beta_pd_nt_pa_z", "beta_fd_pairwise_pa", 
	"beta_fd_pairwise_pa_z", "beta_fd_nt_pa", "beta_fd_nt_pa_z"))
	
bbsgammawithin <- read.csv(file.path(fpdata, 'bbs/biodiversity_CSVs/bbs_withinroute_gamma.csv'), stringsAsFactors = FALSE) %>%
	setDT %>%
	dcast(rteNo ~ radius, value.var = c("richness", "MPD_pa_z", "MNTD_pa_z", "MPDfunc_pa_z", "MNTDfunc_pa_z")) %>%
	setNames(c(names(.)[1], paste('gamma', names(.)[-1], sep = '_')))

bbsbiowithin <- bbsalphawithin %>%
	left_join(bbsbetawithin) %>%
	left_join(bbsgammawithin)	

write.csv(bbsbiowithin, file = file.path(fpdata, 'bbs/biodiversity_CSVs/bbs_allbio_withinroute_wide.csv'), row.names = FALSE)
	
#############################################################
# 22 Feb. 2018

# Geodiversity variables in wide format for FIA.
# Must be done in separate data frames because they are so big.

library(dplyr)
library(reshape2)
fpdata <- '/mnt/research/nasabio/data'

# Point values
fiageopt <- read.csv(file.path(fpdata, 'fia/fia_geo_by_point.csv'), stringsAsFactors = FALSE)
fiahuc <- read.csv(file.path(fpdata, 'fia/fia_huc4.csv'), stringsAsFactors = FALSE)

fiageopt <- fiageopt %>%
	left_join(fiahuc) %>%
	select(PLT_CN, HUC4, everything()) %>%
	setNames(c(names(.)[1:2], paste(names(.)[-(1:2)], 'point', sep = '_')))
	
write.csv(fiageopt, file.path(fpdata, 'fia/fia_geo_by_point.csv'), row.names = FALSE) # Overwrite with HUC4 added.

# Radius values
var_categ <- c('bio1k','bio5k','biocloud','elev','other')

csv_names <- file.path(fpdata, 'fia/geodiv', paste0('fia_usa_', var_categ, '.csv'))

read_reduced <- function(x) {
	read.csv(x, stringsAsFactors = FALSE) %>% select(PLT_CN, radius, mean, sd, variable, richness_geodiv, diversity_geodiv, mode) 
}

fiageorad <- lapply(csv_names, read_reduced)
fiageorad <- lapply(fiageorad, function(x) setNames(x, gsub('nighlight', 'nightlight', names(x)))) # correct typo	

discrete_vars <- grep('soil|geol', unique(fiageorad[[5]]$variable), value = TRUE)

get_stat <- function(x, stat) {
	x %>%
		dcast(PLT_CN ~ variable + radius, value.var = stat) %>%
		setNames(c(names(.)[1], paste(names(.)[-1], stat, sep = '_')))
}

fiageorad_discrete <- filter(fiageorad[[5]], variable %in% discrete_vars)
fiageorad[[5]] <- filter(fiageorad[[5]], !variable %in% discrete_vars)

fiageorad_wide_sd <- lapply(fiageorad, function(x) {
	tri_vars <- grep('_tri|_roughness', unique(x$variable), value = TRUE)
	get_stat(x = filter(x, !variable %in% tri_vars), stat = 'sd')
})
fiageorad_wide_mean <- lapply(fiageorad, get_stat, stat = 'mean')
fiageorad_wide_div <- get_stat(fiageorad_discrete, stat = 'diversity_geodiv')
fiageorad_wide_rich <- get_stat(fiageorad_discrete, stat = 'richness_geodiv')
fiageorad_wide_mode <- get_stat(fiageorad_discrete, stat = 'mode')

fiageorad_wide_sd[[5]] <- fiageorad_wide_sd[[5]] %>% left_join(fiageorad_wide_div) %>% left_join(fiageorad_wide_rich)
fiageorad_wide_mean[[5]] <- fiageorad_wide_mean[[5]] %>% left_join(fiageorad_wide_mode)

for (i in 1:length(var_categ)) {

	write.csv(fiageorad_wide_sd[[i]], file.path(fpdata, 'fia/geodiv' , paste0('fia_', var_categ[i], '_sd_wide.csv')), row.names = FALSE)
	write.csv(fiageorad_wide_mean[[i]], file.path(fpdata, 'fia/geodiv' , paste0('fia_', var_categ[i], '_mean_wide.csv')), row.names = FALSE)
}

# FIA geodiversity 1 km edition

fiageo1km <- read.csv(file.path(fpdata, 'fia/geodiversity_CSVs/fia_geodiversity_1kmradius.csv'), stringsAsFactors = FALSE)

fiageo1km <- fiageo1km %>%
	select(PLT_CN, radius, mean, sd, variable, richness_geodiv, diversity_geodiv, mode) 

discrete_vars <- grep('soil|geol', unique(fiageo1km$variable), value = TRUE)
tri_vars <- grep('_tri|_roughness', unique(fiageo1km$variable), value = TRUE)
	
# Create wide format.
fiageo1km_wide_sd <- fiageo1km %>% filter(!variable %in% c(discrete_vars, tri_vars)) %>% dcast(PLT_CN ~ variable + radius, value.var = 'sd') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'sd', sep = '_')))
fiageo1km_wide_mean <- fiageo1km %>% filter(!variable %in% discrete_vars) %>% dcast(PLT_CN ~ variable + radius, value.var = 'mean') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mean', sep = '_')))
fiageo1km_wide_div <- fiageo1km %>% filter(variable %in% discrete_vars) %>% dcast(PLT_CN ~ variable + radius, value.var = 'diversity_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'diversity', sep = '_')))
fiageo1km_wide_rich <- fiageo1km %>% filter(variable %in% discrete_vars) %>% dcast(PLT_CN ~ variable + radius, value.var = 'richness_geodiv') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'richness', sep = '_')))
fiageo1km_wide_mode <- fiageo1km %>% filter(variable %in% discrete_vars) %>% dcast(PLT_CN ~ variable + radius, value.var = 'mode') %>% setNames(c(names(.)[1], paste(names(.)[-1], 'mode', sep = '_')))

# Join all

fiadat <- fiageo1km_wide_sd %>%
	left_join(fiageo1km_wide_mean) %>%
	left_join(fiageo1km_wide_div) %>%
	left_join(fiageo1km_wide_rich) %>%
	left_join(fiageo1km_wide_mode) %>%
  setNames(gsub('_rtri_', '_tri_', names(.))) # correct typo

	
write.csv(fiadat, file = file.path(fpdata, 'fia/geodiversity_CSVs/fia_allgeo_wide_1kmradius.csv'), row.names = FALSE)

###########################################

# FIA biodiversity

fiaalphapt <- read.csv(file.path(fpdata, 'fia/fiausa_alphadiv.csv'), stringsAsFactors = FALSE) %>% 
	setNames(c(names(.)[1], paste('alpha', names(.)[-(1)], 'point', sep = '_')))

library(data.table)

# Radius values
fiaalpharadius <- read.csv(file.path(fpdata, 'fia/fiausa_alpha.csv'), stringsAsFactors = FALSE)
fiabetaradius <- read.csv(file.path(fpdata, 'fia/fiausa_betatd.csv'), stringsAsFactors = FALSE)
fiagammaradius <- read.csv(file.path(fpdata, 'fia/fiausa_gamma.csv'), stringsAsFactors = FALSE)
n_alpha <- names(fiaalpharadius)[-(1:2)]
n_beta <- names(fiabetaradius)[-(1:2)]
n_gamma <- names(fiagammaradius)[-(1:2)]

fiaalpharadius <- fiaalpharadius %>%
	setDT %>%
	dcast(PLT_CN ~ radius, value.var = n_alpha) %>%
	setNames(c(names(.)[1], paste('alpha', names(.)[-1], sep = '_')))

fiabetaradius <- fiabetaradius %>%	
	setDT %>%
	dcast(PLT_CN ~ radius, value.var = n_beta)
	
fiagammaradius <- fiagammaradius %>%
	setDT %>%
	dcast(PLT_CN ~ radius, value.var = n_gamma) %>%
	setNames(c(names(.)[1], paste('gamma', names(.)[-1], sep = '_')))
	
fiabiodat <- fiaalphapt %>%
	left_join(fiaalpharadius) %>%
	left_join(fiabetaradius) %>%
	left_join(fiagammaradius)
	
write.csv(fiabiodat, file = file.path(fpdata, 'fia/fia_allbio_wide.csv'), row.names = FALSE)