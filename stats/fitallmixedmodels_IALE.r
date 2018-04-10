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

# Edited 08 Mar 2018: Separate the mapping script from the model fitting script.
# New version created for IALE talk on 05 Apr 2018: only do HUC as region, only use the 3 variables from IALE
# However, we will need to do the different radii.
# Edited on 08 Apr 2018: Do linear models with elevation only (3 variables are too complicated)

# Define functions --------------------------------------------------------

# See this old thread for how I composed the model formula: http://r.789695.n4.nabble.com/lmer-with-random-slopes-for-2-or-more-first-level-factors-td902592.html
# Formula should look like:
# lmer(DV ~ IV1 + IV2 + (1|Subject) + (IV1 - 1| Subject) + (IV2 - 1| Subject)) 
# This assumes random slope of IV1 isn't correlated with the random slope of IV2

fit_mm <- function(pred_df, resp_df, pred_vars, resp_var, id_var, region_var, distribution = c('normal', 'beta')) {
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

# Load biodiversity and geodiversity data ---------------------------------

library(dplyr)
fp <- '/mnt/research/nasabio/data'

# BBS
bbsbio <- read.csv(file.path(fp, 'bbs/bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs/bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

# Use 50,100,200 km radius.
bbsbio <- bbsbio %>% 
  select(rteNo, contains('alpha_richness'), contains('beta_td_sorensen_pa'), contains('gamma_richness')) %>%
  select(rteNo, contains('_50'), contains('_100'), contains('_200'), -contains('_500'))
           
                          
bbssp <- bbsgeo[,c('rteNo','lat','lon')]
bbsgeo <- bbsgeo %>%
  select(-contains('richness'), -contains('night')) %>%
  select(rteNo, HUC4, 
         contains('bio12_'), 
         contains('elevation'), 
         contains('footprint')) %>%
  select(rteNo, HUC4, contains('5k')) %>%
  select(rteNo, HUC4, contains('sd'), contains('roughness'), contains('tri')) %>%
  select(rteNo, HUC4, contains('_50_'), contains('_100_'), contains('_200_')) %>%
  mutate(HUC4 = if_else(nchar(HUC4) == 3, paste0('0',HUC4), as.character(HUC4)))

# FIA (only needed variables are in these dfs, rest are on server)
# Use 5,10,20,50,100,200 radius
fiabio <- read.csv(file.path(fp, 'fia/fia_bio_formixedmodels_IALE.csv'), stringsAsFactors = FALSE) 
fiageo <- read.csv(file.path(fp, 'fia/fia_geo_formixedmodels_IALE.csv'), stringsAsFactors = FALSE)
fiageo <- fiageo %>%
  mutate(HUC4 = if_else(nchar(HUC4) == 3, paste0('0',HUC4), as.character(HUC4)))

# Get rid of political border plots and thin down plots -----------------------------------

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
# Also get rid of the plantation plots.
fiasp <- read.csv('~/data/allfia.csv')
plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv', stringsAsFactors = FALSE)
natural <- plantation$PLT_CN[!plantation$plantation]
fiasp <- subset(fiasp, CN %in% fiabio$PLT_CN & CN %in% natural)
fia_aea <- SpatialPoints(coords=fiasp[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
fia_aea <- spTransform(fia_aea, CRS(aea_crs))

fiacoast <- flag_coast_plots(fia_aea, radius = 50e3, border_countries = c('Canada','Mexico')) # Takes ~10 minutes. Try not to repeat this.

fia_aea_noedge <- fia_aea[!fiacoast$is_edge]
fiasp_noedge <- fiasp[!fiacoast$is_edge, ]
# Thin FIA plots to minimum 10 km separation
set.seed(38182)
fiasub10k <- SRS_iterative_N1(fia_aea_noedge, radius = 10e3, n = 20e3, point = sample(length(fia_aea_noedge),1), show_progress = TRUE)
fiaspsub <- fiasp_noedge[fiasub10k, ]
fiabio <- subset(fiabio, PLT_CN %in% fiaspsub$CN)
fiageo <- subset(fiageo, PLT_CN %in% fiaspsub$CN)

# Unfortunately, get rid of the few plots where beta is 0 or 1
fiabio <- fiabio %>%
  mutate_at(vars(contains('beta_td')), function(x) if_else(x > 0 & x < 1, x, as.numeric(NA)))

# Convert all to long format.
library(reshape2)
library(tidyr)
bbsbio_long <- melt(bbsbio, id.vars = 'rteNo') %>%
  separate(variable, into = c('variable', 'radius'), sep = '_(?=[^_]+$)') # Split on last _

bbsgeo_long <- melt(bbsgeo, id.vars = c('rteNo','HUC4')) %>%
  mutate(variable = gsub('human_footprint', 'humanfootprint', variable)) %>%
  separate(variable, into = c('variable', 'metric'), extra = 'merge') %>%
  separate(metric, into = c('metric', 'summary_stat'), sep = '_(?=[^_]+$)') %>%
  separate(metric, into = c('metric', 'radius'), sep = '_(?=[^_]+$)') %>%
  separate(metric, into = c('resolution','metric'), fill = 'right') %>%
  mutate(metric = if_else(is.na(metric), 'sd', metric))

fiabio_long <- melt(fiabio, id.vars = c('PLT_CN','radius'))

fiageo_long <- melt(fiageo, id.vars = c('PLT_CN','HUC4')) %>%
  mutate(variable = gsub('human_footprint', 'humanfootprint', variable)) %>%
  separate(variable, into = c('variable', 'metric'), extra = 'merge') %>%
  separate(metric, into = c('metric', 'summary_stat'), sep = '_(?=[^_]+$)') %>%
  separate(metric, into = c('metric', 'radius'), sep = '_(?=[^_]+$)') %>%
  separate(metric, into = c('resolution','metric'), fill = 'right') %>%
  mutate(metric = if_else(is.na(metric), 'sd', metric))

save(bbsbio_long, bbsgeo_long, fiabio_long, fiageo_long, file = '/mnt/research/nasabio/temp/mm_iale_data.RData')

# Fit models to BBS and FIA data ------------------------------------------

# Factors:
# Metric (SD, roughness, TRI) 3
# Scale (radius) 3 for bbs, 6 for fia
# Diversity type 3

# 27 BBS models, 54 FIA models

load('/mnt/research/nasabio/temp/mm_iale_data.RData')

# Fit all BBS response variables to same predictors by HUC
# By scale and response variable. Fit all predictors and do model selection

library(purrr)
library(dplyr)
library(reshape2)

# Create predictor matrix for each radius
bbs_predictors <- bbsgeo_long %>%
  mutate(metric = if_else(metric=='sd', 'raw', metric),
         radius = as.numeric(radius)) %>%
  group_by(radius) %>%
  do(preds = dcast(., rteNo + HUC4 ~ variable + metric + summary_stat) %>%
       select('rteNo', 'HUC4', contains('elevation')))
  
fia_predictors <- fiageo_long %>%
  mutate(metric = if_else(metric=='sd', 'raw', metric),
         radius = as.numeric(radius)) %>%
  group_by(radius) %>%
  do(preds = dcast(., PLT_CN + HUC4 ~ variable + metric + summary_stat) %>%
       select('PLT_CN', 'HUC4', contains('elevation')))
  
# Create response variables for each radius
bbs_response <- bbsbio_long %>%
  mutate(radius = as.numeric(radius)) %>%
  group_by(radius) %>%
  do(resps = dcast(., rteNo ~ variable))

fia_response <- fiabio_long %>%
  mutate(radius = as.numeric(radius)) %>%
  group_by(radius) %>%
  do(resps = dcast(., PLT_CN ~ variable))


# Do a very basic thing without the mixed model.
# Just to capture the pattern.

pred_resp_combo <- expand.grid(response = c('alpha','beta','gamma'), 
                               predictor = c('sd', 'roughness', 'tri'))

library(betareg)

fit_lm_coefs <- function(dat, resp, pred) {
  dat <- dat[, c(grep(resp, names(dat)), grep(pred, names(dat)))]
  # Make sure predictors are scaled
  dat[,-1] <- scale(dat[,-1])
  model_formula <- formula(paste(names(dat)[1], '~', paste(names(dat)[-1], collapse = '+')))
  if (resp == 'beta') {
    fit <- betareg(model_formula, dat, link = 'logit')
    r2 <- summary(fit)$pseudo.r.squared
    boot_fn <- function(boot_dat, boot_idx) {
      fit <- betareg(model_formula, boot_dat[boot_idx, ], link = 'logit')
      coef(fit)[2]
    }
  } else {
    fit <- lm(model_formula, dat)
    r2 <- summary(fit)$r.squared
    boot_fn <- function(boot_dat, boot_idx) {
      fit <- lm(model_formula, boot_dat[boot_idx, ])
      coef(fit)[2]
    }
  }
  coefs <- coef(fit)
  
  # Include bootstrap.
  require(boot)
  boot_coef <- boot(dat, boot_fn, R = 99, sim = 'ordinary')
  
  #data.frame(slope_precip = coefs[2], slope_elev = coefs[3], slope_footprint = coefs[4], r2 = r2)
  data.frame(slope_elev = coefs[2], r2 = r2, slope_stderr = sd(boot_coef$t))
  
}

coef_by_combo <- function(x) {
  pred_resp_combo %>%
    rowwise %>%
    do(cbind(data.frame(response = .$response, predictor = .$predictor, fit_lm_coefs(x, resp = .$response, pred = .$predictor))))
}

bbs_dat <- map2(bbs_response$resps, bbs_predictors$preds, left_join)
bbs_coefs <- map(bbs_dat, coef_by_combo)
bbs_coefs <- cbind(radius = rep(c(50,100,200), each=9), do.call(rbind, bbs_coefs))

fia_dat <- map2(fia_response$resps, fia_predictors$preds, left_join)
fia_coefs <- map(fia_dat, coef_by_combo)
fia_coefs <- cbind(radius = rep(c(5,10,20,50,100,200), each=9), do.call(rbind, fia_coefs))

write.csv(bbs_coefs, '/mnt/research/nasabio/temp/bbs_coefs_iale.csv', row.names = FALSE)
write.csv(fia_coefs, '/mnt/research/nasabio/temp/fia_coefs_iale.csv', row.names = FALSE)