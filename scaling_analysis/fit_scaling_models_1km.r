# Fit models for scaling analysis: version 2
# QDR NASABioXGeo 1 July 2018

# Modified 15 August 2018: try to force the intercept through zero to aid convergence.

# 1 km response variable, alpha diversity only, using mixed scale sets of predictors

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

args <- commandArgs(TRUE)

for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))
}

NC <- 3
# If arguments are not specified, give default values
if (!exists('NI')) NI <- 5000
if (!exists('NW')) NW <- 3000
if (!exists('delta')) delta <- 0.9

# (4 predictor groups * 3 scales + 1 null) = 13 models * 2 taxa * 2 types of random effect (spatial vs nonspatial) * 2 response variables = 104 models

# Assign options for the model fit being done in this task
climate_prednames <- c('bio1_1k_', 'bio1_1k_tri_', 'bio12_1k_', 'bio12_1k_tri_')
geo_prednames <- c('soil_type_5k_', 'geological_age_1k_', 'dhi_gpp_1k_', 'dhi_gpp_1k_tri_', 'elevation_30m_', 'elevation_30m_tri_')
human_prednames <- c('human_footprint_1k_', 'human_footprint_1k_tri_')

# Mixed scale predictor sets.
# 250 combinations, including all scales
task_table <- expand.grid(taxon = c('fia','bbs'),
                         rv = c('alpha_richness'), # Do only alpha
                         random_effect = c('spatial'), # Do only the spatial random effect (don't care about testing the nonspatial random effect)
						 climate_scale = c(0, 1, 5, 20, 100),
						 geo_scale = c(0, 1, 5, 20, 100),
						 human_scale = c(0, 1, 5, 20, 100),
                         stringsAsFactors = FALSE)

taxon <- task_table$taxon[task]
climate_scale <- task_table$climate_scale[task]
geo_scale <- task_table$geo_scale[task]
human_scale <- task_table$human_scale[task]
predictors <- ifelse(climate_scale==0 & geo_scale==0 & human_scale==0, 'null', 'yes')
random_effect <- task_table$random_effect[task]
rv <- task_table$rv[task]

# Function used for model fitting.
source('/mnt/research/nasabio/code/fit_mv_mm.r')

options(mc.cores = 3)

# Load data for appropriate taxon
if (taxon == 'bbs') {
  load('/mnt/research/nasabio/temp/bbs_scaling_dat_1km.RData')
  geodat <- bbsgeo
  biodat <- bbsbio
  siteid <- 'rteNo'


} else {
  load('/mnt/research/nasabio/temp/fia_scaling_dat_1km.RData')
  geodat <- fiageo
  biodat <- fiabio
  siteid <- 'PLT_CN'
  
  # Edit 18 June: Remove one unexplained outlier (~5E14) of elevation TRI
  geodat$elevation_30m_tri_100_mean[geodat$elevation_30m_tri_100_mean > 1e6] <- NA
}

# Get correct subset of predictors by type and scale
prednames <- names(geodat)[!(names(geodat) %in% c('rteNo','HUC4','TNC','BCR'))]

climate_prednames <- grep(paste(paste0(climate_prednames, climate_scale, '_'), collapse='|'), prednames, value = TRUE)
geo_prednames <- grep(paste(paste0(geo_prednames, geo_scale, '_'), collapse='|'), prednames, value = TRUE)
human_prednames <- grep(paste(paste0(human_prednames, human_scale, '_'), collapse='|'), prednames, value = TRUE)

prednames <- c(climate_prednames, geo_prednames, human_prednames)

distrib <- 'gaussian' # Fit all with Gaussian
ecoregion <- 'TNC' # Use TNC as the ecoregion

added_priors <- NULL # For most tasks, model will converge using only default priors

# Edit 18 June: Add priors for the models that did not converge
# Edit 15 Aug: remove this because there is no intercept anymore.
#if (rv == 'alpha_richness' & taxon == 'fia') {
#  added_priors <- c(brms::set_prior('student_t(5, 0, 2)', class = 'Intercept'))
#} 

fit <- fit_mv_mm(pred_df = geodat, 
				  resp_df = biodat, 
				  pred_vars = prednames, 
				  resp_vars = rv, 
				  id_var = siteid, 
				  region_var = ecoregion, 
				  distribution = distrib, 
				  adj_matrix = tnc_bin,
				  priors = added_priors,
				  n_chains = NC,
				  n_iter = NI,
				  n_warmup = NW,
				  delta = delta,
				  random_effect_type = random_effect,
				  force_zero_intercept = TRUE
				  )

# Save all fits
save(fit, file = paste0('/mnt/research/nasabio/temp/scalingfits/1kmresponse_fit',task,'.RData'))
