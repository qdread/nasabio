# Mixed model with random effects for bbs predicting geo by bio
# Right now, just use HUCs as the random effect, ignoring space
# Use random slope and intercept

# Just use a few predictors for alpha, beta, and gamma (modified 07 mar to add human footprint)
# Func, phy, and tax

# Make a map of relationship.

# Use only 100 km radius. Use gamma.
resps <- bbsbio %>% select(rteNo,
                           alpha_richness_100, beta_td_sorensen_pa_100, gamma_richness_100, 
                           alpha_MPD_pa_z_100, beta_pd_pairwise_pa_z_100, gamma_MPD_pa_z_100,
                           alpha_MPDfunc_pa_z_100, beta_fd_pairwise_pa_z_100, gamma_MPDfunc_pa_z_100)

# Each predictor has 3 versions (mean, sd, and tri mean). Point values not used here.
# Soil type and geo age will only have 3 predictors but some are split to dummy vars 
# (8 predictors * 3 ways of aggregating = ~24 predictors . . . .actually 22 + 1 random effect)
# Temp, Precip, Clouds, Elevation, Soil Type, Geo Age, Human Footprint, GPP
preds <- bbsgeo %>% 
  select(rteNo, HUC4, BCR, TNC,
                           contains('bio1_'), contains('bio12_'), contains('biocloud1_'),
                           contains('elevation'), contains('soil'), contains('geological'),
                           contains('footprint'), contains('gpp')) %>%
  select(-contains('point'))



preds_scale <- cbind(preds[,1:4], scale(preds[,!sapply(preds,is.factor)][,-(1:4)]))

gamma_dat <- resps %>%
  select(rteNo, gamma_richness_100) %>%
  left_join(preds_scale)

library(usdm)
vif(preds_scale[,c('elevation_5k_100_mean', 'elevation_5k_100_sd', 'bio1_5k_100_mean' , 'bio1_5k_100_sd' , 'geological_age_5k_100_diversity' , 'soil_type_5k_100_diversity' , 'bio12_5k_100_mean' , 'bio12_5k_100_sd' , 'dhi_gpp_5k_100_mean' , 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean', 'human_footprint_5k_100_sd')])

cor(preds_scale[,c('elevation_5k_100_mean', 'elevation_5k_100_sd', 'bio1_5k_100_mean' , 'bio1_5k_100_sd' , 'geological_age_5k_100_diversity' , 'soil_type_5k_100_diversity' , 'bio12_5k_100_mean' , 'bio12_5k_100_sd' , 'dhi_gpp_5k_100_mean' , 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean', 'human_footprint_5k_100_sd')], use = 'pairwise.complete.obs')

# Mean GPP is highly correlated with precipitation.
# SD of temperature is highly correlated with SD of elevation.
# Mean elevation is highly correlated with SD of elevation.
# Mean human footprint is highly correlated with SD of human footprint
# Removing them makes everything much better.
vif(preds_scale[,c('elevation_5k_100_sd', 'bio1_5k_100_mean'  , 'geological_age_5k_100_diversity' , 'soil_type_5k_100_diversity' , 'bio12_5k_100_mean' , 'bio12_5k_100_sd'  , 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')])

# Do just a few here.
# Linear, additive (no nonlinearities or interactions yet)
library(lme4)

# This one fits random intercepts but no random slopes
mm1 <- lmer(gamma_richness_100 ~ elevation_5k_100_sd + bio1_5k_100_mean + geological_age_5k_100_diversity + soil_type_5k_100_diversity + bio12_5k_100_mean + bio12_5k_100_sd + dhi_gpp_5k_100_sd + (1|HUC4), data = gamma_dat)

# Need to also fit random slopes.
mm2 <- lmer(gamma_richness_100 ~ elevation_5k_100_sd + bio1_5k_100_mean + geological_age_5k_100_diversity + soil_type_5k_100_diversity + bio12_5k_100_mean + bio12_5k_100_sd + dhi_gpp_5k_100_sd + (elevation_5k_100_sd|HUC4) + (bio1_5k_100_mean|HUC4) + (geological_age_5k_100_diversity|HUC4) + (soil_type_5k_100_diversity|HUC4) + (bio12_5k_100_mean|HUC4) + (bio12_5k_100_sd|HUC4) + (dhi_gpp_5k_100_sd|HUC4), data = gamma_dat[complete.cases(gamma_dat), ])

# Just random slopes without fixed slopes
mm3 <- lmer(gamma_richness_100 ~ (bio1_5k_100_mean|HUC4) + (elevation_5k_100_sd|HUC4) + (geological_age_5k_100_diversity|HUC4) + (soil_type_5k_100_diversity|HUC4) + (bio12_5k_100_mean|HUC4) + (bio12_5k_100_sd|HUC4) + (dhi_gpp_5k_100_sd|HUC4), data = gamma_dat[complete.cases(gamma_dat), ])

# Determine how to extract the coefficients here, then make a map
# Also try alternative regionalizations besides HUCs
mm_coef <- coef(mm2)[[1]]
mm_coef$HUC4 <- dimnames(mm_coef)[[1]]
mm_coef$HUC4[nchar(mm_coef$HUC4) == 3] <- paste0('0', mm_coef$HUC4[nchar(mm_coef$HUC4) == 3])

mm2_conf <- confint(mm2, method = 'Wald') # Gives confidence intervals for fixed effects only.

# Load HUC4 boundary map
library(sp)
library(rgdal)
library(ggplot2)
huc4 <- readOGR(dsn = 'C:/Users/Q/Dropbox/projects/aquaxterra/hucshapefiles', layer = 'HU4_CONUS_Alb')
huc4@data <- huc4@data %>% 
  mutate(id = rownames(huc4@data), HUC4 = as.character(HUC4)) %>% 
  left_join(mm_coef)

huc4_fort <- fortify(huc4, region = 'id') %>% left_join(huc4@data, by = 'id')

# Overall color scale for all
coefrange <- sapply(mm_coef[,2:8], range)
coefmax <- ceiling(max(abs(coefrange)))

rbcolors <- rev(RColorBrewer::brewer.pal(9, 'RdYlBu'))

vars_to_plot <- names(mm_coef)[2:8]

states <- read.csv('~/R/states_albers.csv', stringsAsFactors = FALSE)

for (i in vars_to_plot) {
  map_i <- ggplot(huc4_fort) +
    geom_polygon(aes_string(x='long', y='lat', group='group', fill=i)) +
    geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
    geom_path(data = states, aes(x = long, y = lat, group = group), color = 'gray20') +
    scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0) +
    coord_equal() +
    theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'black'), panel.border = element_blank(), plot.background = element_rect(fill = 'black'), legend.position = c(0.13,0.1), legend.direction = 'horizontal', legend.title = element_blank())
  ggsave(filename = paste0('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_coefficient_maps/huc4_', i, '_map.png'), plot = map_i, height = 6, width = 9, dpi = 400)
  print(i)
}


# Alternative regions:

# Bird conservation regions

# Fit model with BCR
mm_bcr <- lmer(gamma_richness_100 ~ bio1_5k_100_mean + elevation_5k_100_sd + geological_age_5k_100_diversity + soil_type_5k_100_diversity + bio12_5k_100_mean + bio12_5k_100_sd + dhi_gpp_5k_100_sd + (bio1_5k_100_mean|BCR) + (elevation_5k_100_sd|BCR) + (geological_age_5k_100_diversity|BCR) + (soil_type_5k_100_diversity|BCR) + (bio12_5k_100_mean|BCR) + (bio12_5k_100_sd|BCR) + (dhi_gpp_5k_100_sd|BCR), data = gamma_dat[complete.cases(gamma_dat), ])

coef_bcr <- coef(mm_bcr)[[1]]
coef_bcr$BCR <- dimnames(coef_bcr)[[1]]

# Make maps with BCR
bcr <- readOGR(dsn = 'C:/Users/Q/Dropbox/projects/nasabiodiv/regions', layer = 'BCR_Terrestrial_master')
# Reproject bcr shapefile and get rid of the ones outside the continental USA.
bcr <- subset(bcr, BCRNAME %in% coef_bcr$BCR)
bcr@data <- bcr@data %>% 
  mutate(id = rownames(bcr@data), BCR = as.character(BCRNAME)) %>% 
  left_join(coef_bcr)

bcr <- spTransform(bcr, CRSobj = CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))


bcr_fort <- fortify(bcr, region = 'id') %>% left_join(bcr@data, by = 'id')

rbcolors <- rev(RColorBrewer::brewer.pal(9, 'RdYlBu'))

vars_to_plot <- names(coef_bcr)[2:8]

states <- read.csv('~/R/states_albers.csv', stringsAsFactors = FALSE)

for (i in vars_to_plot) {
  map_i <- ggplot(bcr_fort) +
    geom_polygon(aes_string(x='long', y='lat', group='group', fill=i)) +
    geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
    geom_path(data = states, aes(x = long, y = lat, group = group), color = 'gray20') +
    scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0) +
    coord_equal() +
    theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'black'), panel.border = element_blank(), plot.background = element_rect(fill = 'black'), legend.position = c(0.13,0.1), legend.direction = 'horizontal', legend.title = element_blank())
  ggsave(filename = paste0('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_coefficient_maps/BCR_', i, '_map.png'), plot = map_i, height = 6, width = 9, dpi = 400)
  print(i)
}

# TNC ecoregions

# Fit model with TNC
mm_tnc <- lmer(gamma_richness_100 ~ bio1_5k_100_mean + elevation_5k_100_sd + geological_age_5k_100_diversity + soil_type_5k_100_diversity + bio12_5k_100_mean + bio12_5k_100_sd + dhi_gpp_5k_100_sd + (bio1_5k_100_mean|TNC) + (elevation_5k_100_sd|TNC) + (geological_age_5k_100_diversity|TNC) + (soil_type_5k_100_diversity|TNC) + (bio12_5k_100_mean|TNC) + (bio12_5k_100_sd|TNC) + (dhi_gpp_5k_100_sd|TNC), data = gamma_dat[complete.cases(gamma_dat), ])

coef_tnc <- coef(mm_tnc)[[1]]
coef_tnc$TNC <- dimnames(coef_tnc)[[1]]

# Make maps with TNC
tnc <- readOGR(dsn = 'C:/Users/Q/Dropbox/projects/nasabiodiv/regions', layer = 'tnc_terr_ecoregions')
# Reproject tnc shapefile and get rid of the ones outside the continental USA.
tnc <- subset(tnc, ECODE_NAME %in% coef_tnc$TNC)
tnc@data <- tnc@data %>% 
  mutate(id = rownames(tnc@data), TNC = as.character(ECODE_NAME)) %>% 
  left_join(coef_tnc)

tnc <- spTransform(tnc, CRSobj = CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

tnc_fort <- fortify(tnc, region = 'id') %>% left_join(tnc@data, by = 'id')

rbcolors <- rev(RColorBrewer::brewer.pal(9, 'RdYlBu'))

vars_to_plot <- names(coef_tnc)[2:8]

states <- read.csv('~/R/states_albers.csv', stringsAsFactors = FALSE)

for (i in vars_to_plot) {
  map_i <- ggplot(tnc_fort) +
    geom_polygon(aes_string(x='long', y='lat', group='group', fill=i)) +
    geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
    geom_path(data = states, aes(x = long, y = lat, group = group), color = 'gray20') +
    scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0) +
    coord_equal() +
    theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'black'), panel.border = element_blank(), plot.background = element_rect(fill = 'black'), legend.position = c(0.13,0.1), legend.direction = 'horizontal', legend.title = element_blank())
  ggsave(filename = paste0('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_coefficient_maps/TNC_', i, '_map.png'), plot = map_i, height = 6, width = 9, dpi = 400)
  print(i)
}