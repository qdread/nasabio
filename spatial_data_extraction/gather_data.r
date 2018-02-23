
# 21 Feb 2018
# All geodiversity variables in a single wide format data frame for BBS and FIA.

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
	left_join(bbsgeorad_wide_mode)

# correct typo
bbsdat <- setNames(bbsdat, gsub('nighlight', 'nightlight', names(bbsdat)))
	
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