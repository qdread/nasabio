# Mixed model with non-spatial random effects
# General function that takes the following arguments:
# Dataset (bbs or fia)
# Response variable (a/b/g and t/f/p diversity)
# Predictor variables (narrowed down ahead of time)
# Region variable to use as random effect (huc4, bcr, or tnc)
# ID variable is rteNo for BBS and PLT_CN for FIA

# One function to fit model and return coefficients
# Another function to create and save the map

# QDR NASABIOXGEO 07 Mar 2018

# Define functions --------------------------------------------------------

# See this old thread for how I composed the model formula: http://r.789695.n4.nabble.com/lmer-with-random-slopes-for-2-or-more-first-level-factors-td902592.html
# Formula should look like:
# lmer(DV ~ IV1 + IV2 + (1|Subject) + (IV1 - 1| Subject) + (IV2 - 1| Subject)) 
# This assumes random slope of IV1 isn't correlated with the random slope of IV2

fit_mm <- function(pred_df, resp_df, pred_vars, resp_var, id_var, region_var, distribution =c('normal','beta')) {
  require(lme4)
  require(dplyr)
  require(glmmTMB)
  # Build formula and data
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_var)]
  pred_df <- pred_df[, c(id_var, pred_vars)]
  pred_df[,-1] <- scale(pred_df[,-1])
  pred_var_names <- names(pred_df)[-1]
  region_name <- names(id_df)[2]
  fixed_effects <- paste(pred_var_names, collapse = '+')
  random_effects <- paste(c(paste('(1|', region_var, ')', sep = ''), paste('(', pred_var_names, ' - 1|', region_var, ')', sep = '')), collapse = '+')
  formula_string <- paste(names(resp_df)[2], '~', fixed_effects, '+', random_effects)
  dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  # Fit model, extract coefficients, and format them
  # Use lmer() if normally distributed response, use glmmTMB() if beta distributed response
  if (distribution == 'normal') {
    mm <- lmer(formula(formula_string), data = dat)
    mm_coef <- coef(mm)[[1]]
    mm_coef <- cbind(region = dimnames(mm_coef)[[1]], mm_coef)
    mm_coef <- mm_coef[,-grep('Intercept', names(mm_coef))]
  }
  if (distribution == 'beta') {
    mm <- glmmTMB(formula(formula_string), data = dat, family = list(family = 'beta', link = 'logit'))
    mm_coef <- ranef(mm)[[1]][[1]]
    mm_coef <- cbind(region = dimnames(mm_coef)[[1]], mm_coef)
    mm_coef <- mm_coef[,-grep('Intercept', names(mm_coef))]
  }
  return(list(model = mm, coef = mm_coef))
}

# Make map (run fortify each time this function is called, slower but simpler)
# Return the plots as a list so that we can tile them in different ways
model_map <- function(coefs, vars_to_plot, fill_scale, regions, state_borders) {
  regions@data <- regions@data %>% 
    left_join(coefs)
  
  region_fort <- fortify(regions, region = 'id') %>% left_join(regions@data, by = 'id')
  
  blktheme <- theme_bw() + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          panel.grid = element_blank(), 
          panel.background = element_rect(color = 'black', fill = 'black'), 
          panel.border = element_blank(), 
          plot.background = element_rect(fill = 'black'), 
          legend.position = c(0.13,0.1), 
          legend.direction = 'horizontal', 
          legend.title = element_blank())
  
  map_list <- list()
  for (i in vars_to_plot) {
    map_list[[length(map_list) + 1]] <- ggplot(region_fort) +
      geom_polygon(aes_string(x='long', y='lat', group='group', fill=i)) +
      geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
      geom_path(data = state_borders, aes(x = long, y = lat, group = group), color = 'gray20') +
      fill_scale +
      coord_equal() +
      blktheme
    print(i)
  }
  return(map_list)
}

# Make predicted vs observed plot
obs_pred_plot <- function() {
  
}

# Load biodiversity and geodiversity data ---------------------------------

library(dplyr)
fp <- '/mnt/research/nasabio/data'

# BBS
bbsbio <- read.csv(file.path(fp, 'bbs/bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs/bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

# Use only 100 km radius.
bbsbio <- bbsbio %>% select(rteNo,
                           alpha_richness_100, beta_td_sorensen_pa_100, gamma_richness_100, 
                           alpha_MPD_pa_z_100, beta_pd_pairwise_pa_z_100, gamma_MPD_pa_z_100,
                           alpha_MPDfunc_pa_z_100, beta_fd_pairwise_pa_z_100, gamma_MPDfunc_pa_z_100)

bbsgeo <- bbsgeo %>%
  select(rteNo, HUC4, BCR, TNC, contains('point'), matches('5k.*100')) %>%
  select(-contains('roughness'), -contains('richness'), -contains('night')) %>%
  select(rteNo, HUC4, BCR, TNC,
         contains('bio1_'), contains('bio12_'), contains('biocloud1_'),
         contains('elevation'), contains('soil'), contains('geological'),
         contains('footprint'), contains('gpp')) %>%
  select(-contains('point')) %>%
  mutate(HUC4 = if_else(nchar(HUC4) == 3, paste0('0',HUC4), as.character(HUC4)))

# FIA (only needed variables are in these dfs, rest are on server)
# Use effective species number for Shannon
fiabio <- read.csv(file.path(fp, 'fia/fia_bio_formixedmodels.csv'), stringsAsFactors = FALSE) %>%
  mutate(alpha_shannon = exp(alpha_shannon), gamma_shannon = exp(gamma_shannon)) %>%
  rename(alpha_effspn = alpha_shannon, gamma_effspn = gamma_shannon)
fiageo <- read.csv(file.path(fp, 'fia/fia_geo_formixedmodels.csv'), stringsAsFactors = FALSE)
fiageo <- fiageo %>%
  setNames(gsub('_geodiv','',names(.))) %>%
  mutate(HUC4 = if_else(nchar(HUC4) == 3, paste0('0',HUC4), as.character(HUC4)))

# Load region data --------------------------------------------------------

library(sp)
library(rgdal)
library(ggplot2)
fpregion <- file.path(fp, 'ecoregions')
huc4 <- readOGR(dsn = fpregion, layer = 'HU4_CONUS_Alb')
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
states <- read.csv('~/states_albers.csv', stringsAsFactors = FALSE)

# Convert all to Albers with "region" in the name.
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
bcr <- spTransform(bcr, CRSobj = CRS(aea_crs))
tnc <- spTransform(tnc, CRSobj = CRS(aea_crs))

huc4@data <- huc4@data %>%
  mutate(id = rownames(huc4@data), region = as.character(HUC4))

bcr@data <- bcr@data %>%
  mutate(id = rownames(bcr@data), region = as.character(BCRNAME))

tnc@data <- tnc@data %>%
  mutate(id = rownames(tnc@data), region = as.character(ECODE_NAME))
  

# Fit models to BBS and FIA data ------------------------------------------

prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')

# Fit all BBS response variables to same predictors by HUC, BCR, TNC

library(purrr)

distribs <- ifelse(grepl('beta_td', names(bbsbio)[-1]), 'beta', 'normal') # Beta_td follows beta distribution, rest are normally distributed

fit_bbs_huc <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'HUC4', distr))
fit_bbs_bcr <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'BCR', distr))
fit_bbs_tnc <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'TNC', distr))


# Fit all FIA response variables to same predictors by the 3 regions also

distribs <- ifelse(grepl('beta_td', names(fiabio)[-1]), 'beta', 'normal')

fit_fia_huc <- map2(names(fiabio)[-1], distribs, function(varname, distr) fit_mm(fiageo, fiabio, prednames, varname, 'PLT_CN', 'HUC4', distr))
fit_fia_bcr <- map2(names(fiabio)[-1], distribs, function(varname, distr) fit_mm(fiageo, fiabio, prednames, varname, 'PLT_CN', 'BCR', distr))
fit_fia_tnc <- map2(names(fiabio)[-1], distribs, function(varname, distr) fit_mm(fiageo, fiabio, prednames, varname, 'PLT_CN', 'TNC', distr))


# Create all maps ------------------------------------------

rbfill <- scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0)

# Subset out the regions that are outside the US.
bcr <- subset(bcr, region %in% fit_bbs_bcr[[1]][[2]]$region & COUNTRY %in% 'USA' & !PROVINCE_S %in% 'ALASKA')
tnc <- subset(tnc, region %in% fit_bbs_tnc[[1]][[2]]$region)

# Create nested list of maps for each predictor by response combo (9 predictors x 8 responses)
maps_bbs_huc <- map(fit_bbs_huc, function(fit) model_map(fit$coef, prednames, rbfill, huc4, states))
maps_bbs_bcr <- map(fit_bbs_bcr, function(fit) model_map(fit$coef, prednames, rbfill, bcr, states))
maps_bbs_tnc <- map(fit_bbs_tnc, function(fit) model_map(fit$coef, prednames, rbfill, tnc, states))

names(maps_bbs_huc) <- names(bbsbio)[-1]
maps_bbs_huc <- lapply(maps_bbs_huc, setNames, nm = prednames)
names(maps_bbs_bcr) <- names(bbsbio)[-1]
maps_bbs_bcr <- lapply(maps_bbs_bcr, setNames, nm = prednames)
names(maps_bbs_tnc) <- names(bbsbio)[-1]
maps_bbs_tnc <- lapply(maps_bbs_tnc, setNames, nm = prednames)

# Create grids of maps by predictor, and grids of maps by response
library(gridExtra)
fpfig <- '/mnt/research/nasabio/figs/bbs_coefficient_maps'
bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
geo_names <- c('elevation_sd','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','precip_sd','gpp_sd','footprint_mean')
tw <- theme(plot.title = element_text(color = 'white'))
for (i in 1:8) {
  print(i)
  png(file.path(fpfig, paste0('HUC4_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_bbs_huc, bio_titles, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), nrow = 3)
  dev.off()
  png(file.path(fpfig, paste0('BCR_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_bbs_bcr, bio_titles, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), nrow = 3)
  dev.off()
  png(file.path(fpfig, paste0('TNC_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_bbs_tnc, bio_titles, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), nrow = 3)
  dev.off()
}
