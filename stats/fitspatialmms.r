# Mixed model with non-spatial random effects
# General function that takes the following arguments:
# Dataset (bbs or fia)
# Response variable (a/b/g and t/f/p diversity)
# Predictor variables (narrowed down ahead of time)
# Region variable to use as random effect (huc4, bcr, or tnc)
# Adjacency matrix is needed to specify the CAR autocorrelation structure for the random effects.
# ID variable is rteNo for BBS and PLT_CN for FIA

# QDR NASABIOXGEO 24 Apr 2018

# New version of this script created 24 Apr: use brms to fit the spatial model.

# Define functions --------------------------------------------------------

# See this old thread for how I composed the model formula: http://r.789695.n4.nabble.com/lmer-with-random-slopes-for-2-or-more-first-level-factors-td902592.html
# Formula should look like:
# lmer(DV ~ IV1 + IV2 + (1|Subject) + (IV1 - 1| Subject) + (IV2 - 1| Subject)) 
# This assumes random slope of IV1 isn't correlated with the random slope of IV2

fit_spatial_mm <- function(pred_df, resp_df, pred_vars, resp_var, id_var, region_var, adj_matrix, distribution =c('gaussian','beta'), n_chains = 2, n_iter = 2000, n_warmup = 1000) {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_var)]
  pred_df <- pred_df[, c(id_var, pred_vars)]
  pred_df[,-1] <- scale(pred_df[,-1])
  pred_var_names <- names(pred_df)[-1]
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  fixed_effects <- paste(pred_var_names, collapse = '+')
  random_effects <- paste(c(paste('(1|', region_name, ')', sep = ''), paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
  formula_string <- paste(names(resp_df)[2], '~', fixed_effects, '+', random_effects)
  dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  # Fit model, extract coefficients, and format them
  mm <- brm(formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region, type = 'esicar'),
			 chains = n_chains, iter = n_iter, warmup = n_warmup)
  fixed_effects <- fixef(mm)
  random_effects <- ranef(mm)
  fixed_effects <- melt(fixed_effects, varnames = c('parameter', 'stat'))
  random_effects <- melt(random_effects$region, varnames = c('region', 'parameter', 'stat'))
  mm_coef <- rbind(cbind(effect = 'fixed', region = NA, fixed_effects),
				   cbind(effect = 'random', random_effects))
  return(list(model = mm, coef = mm_coef))
}

# Load biodiversity and geodiversity data ---------------------------------

library(dplyr)
fp <- '/mnt/research/nasabio/data'

# Load adjacency matrices
load(file.path(fp, 'ecoregions/ecoregion_adjacency.RData'))

# BBS
bbsbio <- read.csv(file.path(fp, 'bbs/bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs/bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

# Use only 100 km radius.
bbsbio <- bbsbio %>% select(rteNo,
                           alpha_richness_100, beta_td_sorensen_pa_100, gamma_richness_100, 
                           alpha_MPD_pa_z_100, beta_pd_pairwise_pa_z_100, gamma_MPD_pa_z_100,
                           alpha_MPDfunc_pa_z_100, beta_fd_pairwise_pa_z_100, gamma_MPDfunc_pa_z_100)

bbssp <- bbsgeo[,c('rteNo','lat','lon')]
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

# Get rid of coasts and thin down plots -----------------------------------

library(sp)

aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
bbssp <- subset(bbssp, rteNo %in% bbsbio$rteNo)
bbs_aea <- SpatialPoints(coords=bbssp[,c('lon','lat')], proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
bbs_aea <- spTransform(bbs_aea, CRS(aea_crs))

# Function for iterative search.
source('/mnt/research/nasabio/code/SRS_iterative.r')
# Functions for flagging edge plots
source('/mnt/research/nasabio/code/spatial_fns.r')

# Use 50km to Can/Mex. Don't thin BBS plots
bbscoast <- flag_coast_plots(bbs_aea, radius = 50e3, border_countries = c('Canada','Mexico'))
noedge_rte <- bbssp$rteNo[!bbscoast$is_edge]
bbsbio <- subset(bbsbio, rteNo %in% noedge_rte)
bbsgeo <- subset(bbsgeo, rteNo %in% noedge_rte)

# For FIA, use 10 km thinning and get rid of 50km to Can/Mex. This will reduce dataset to very manageable <10k plots.
fiasp <- read.csv('~/data/allfia.csv')
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
fiasp <- subset(fiasp, CN %in% fiabio$PLT_CN)
fia_aea <- SpatialPoints(coords=fiasp[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
fia_aea <- spTransform(fia_aea, CRS(aea_crs))

fiacoast <- flag_coast_plots(fia_aea, radius = 50e3, border_countries = c('Canada','Mexico'))

fia_aea_noedge <- fia_aea[!fiacoast$is_edge]
fiasp_noedge <- fiasp[!fiacoast$is_edge, ]
# Thin FIA plots to minimum 10 km separation
set.seed(38182)
fiasub10k <- SRS_iterative_N1(fia_aea_noedge, radius = 10e3, n = 2e5, point = sample(length(fia_aea_noedge),1), show_progress = TRUE)
fiaspsub <- fiasp_noedge[fiasub10k, ]
fiabio <- subset(fiabio, PLT_CN %in% fiaspsub$CN)
fiageo <- subset(fiageo, PLT_CN %in% fiaspsub$CN)

# Unfortunately, get rid of the few plots where beta is 0 or 1
fiabio <- fiabio %>%
  mutate_at(vars(contains('beta_td')), function(x) if_else(x > 0 & x < 1, x, as.numeric(NA)))

# Fit models to BBS and FIA data ------------------------------------------

prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')

# Fit all BBS response variables to same predictors by HUC, BCR, TNC

library(purrr)
options(mc.cores = 2)

distribs <- ifelse(grepl('beta_td', names(bbsbio)[-1]), 'beta', 'normal') # Beta_td follows beta distribution, rest are normally distributed

fit_bbs_huc <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'HUC4', distr))
fit_bbs_bcr <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'BCR', distr))
fit_bbs_tnc <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'TNC', distr))


# Fit all FIA response variables to same predictors by the 3 regions also

distribs <- ifelse(grepl('beta_td', names(fiabio)[-1]), 'beta', 'normal')

fit_fia_huc <- map2(names(fiabio)[-1], distribs, function(varname, distr) fit_mm(fiageo, fiabio, prednames, varname, 'PLT_CN', 'HUC4', distr))
fit_fia_bcr <- map2(names(fiabio)[-1], distribs, function(varname, distr) fit_mm(fiageo, fiabio, prednames, varname, 'PLT_CN', 'BCR', distr))
fit_fia_tnc <- map2(names(fiabio)[-1], distribs, function(varname, distr) fit_mm(fiageo, fiabio, prednames, varname, 'PLT_CN', 'TNC', distr))

# Save all fits
save(fit_bbs_huc, fit_bbs_bcr, fit_bbs_tnc, fit_fia_huc, fit_fia_bcr, fit_fia_tnc, file = '/mnt/research/nasabio/temp/mmfits.RData')

