# Mixed model with random effects for bbs predicting geo by bio
# Right now, just use HUCs as the random effect, ignoring space
# Use random slope and intercept

# Just use a few predictors for alpha, beta, and gamma
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
  select(rteNo, HUC4, 
                           contains('bio1_'), contains('bio12_'), contains('biocloud1_'),
                           contains('elevation'), contains('soil'), contains('geological'),
                           contains('footprint'), contains('gpp')) %>%
  select(-contains('point'))



preds_scale <- cbind(preds[,1:2], scale(preds[,!sapply(preds,is.factor)][,-(1:2)]))

gamma_dat <- resps %>%
  select(rteNo, gamma_richness_100) %>%
  left_join(preds_scale)

library(usdm)
vif(preds_scale[,c('bio1_5k_100_mean' , 'bio1_5k_100_sd' , 'geological_age_5k_100_diversity' , 'soil_type_5k_100_diversity' , 'bio12_5k_100_mean' , 'bio12_5k_100_sd' , 'dhi_gpp_5k_100_mean' , 'dhi_gpp_5k_100_sd')])

cor(preds_scale[,c('bio1_5k_100_mean' , 'bio1_5k_100_sd' , 'geological_age_5k_100_diversity' , 'soil_type_5k_100_diversity' , 'bio12_5k_100_mean' , 'bio12_5k_100_sd' , 'dhi_gpp_5k_100_mean' , 'dhi_gpp_5k_100_sd')], use = 'pairwise.complete.obs')

# Mean GPP is highly correlated with precipitation.
# Removing it makes everything much better.
vif(preds_scale[,c('bio1_5k_100_mean' , 'bio1_5k_100_sd' , 'geological_age_5k_100_diversity' , 'soil_type_5k_100_diversity' , 'bio12_5k_100_mean' , 'bio12_5k_100_sd'  , 'dhi_gpp_5k_100_sd')])

# Do just a few here.
# Linear, additive (no nonlinearities or interactions yet)
library(lme4)

# This one fits random intercepts but no random slopes
mm1 <- lmer(gamma_richness_100 ~ bio1_5k_100_mean + bio1_5k_100_sd + geological_age_5k_100_diversity + soil_type_5k_100_diversity + bio12_5k_100_mean + bio12_5k_100_sd + dhi_gpp_5k_100_sd + (1|HUC4), data = gamma_dat)

# Need to also fit random slopes.
mm2 <- lmer(gamma_richness_100 ~ bio1_5k_100_mean + bio1_5k_100_sd + geological_age_5k_100_diversity + soil_type_5k_100_diversity + bio12_5k_100_mean + bio12_5k_100_sd + dhi_gpp_5k_100_sd + (bio1_5k_100_mean|HUC4) + (bio1_5k_100_sd|HUC4) + (geological_age_5k_100_diversity|HUC4) + (soil_type_5k_100_diversity|HUC4) + (bio12_5k_100_mean|HUC4) + (bio12_5k_100_sd|HUC4) + (dhi_gpp_5k_100_sd|HUC4), data = gamma_dat)

# Determine how to extract the coefficients here, then make a map
# Also try alternative regionalizations besides HUCs
mm_coef <- coef(mm2)[[1]]
mm_coef$HUC4 <- dimnames(mm_coef)[[1]]
mm_coef$HUC4[nchar(mm_coef$HUC4) == 3] <- paste0('0', mm_coef$HUC4[nchar(mm_coef$HUC4) == 3])

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
