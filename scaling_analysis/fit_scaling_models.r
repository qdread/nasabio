# Fit models for scaling analysis
# QDR NASABioXGeo 14 June 2018

# Modified 18 June: add priors for the models that did not converge

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
climate_prednames <- c('bio1_', 'bio12_')
geo_prednames <- c('soil_type', 'geological_age', 'dhi_gpp', 'elevation')
human_prednames <- c('human_footprint')

task_table <- expand.grid(taxon = c('fia','bbs'),
                         rv = c('alpha_richness', 'beta_td_sorensen_pa'),
                         random_effect = c('spatial', 'nonspatial'),
						 scale = c(5, 20, 100),
						 predictor_group = c('full', 'climate', 'geo', 'human'),
                         stringsAsFactors = FALSE)
task_table <- merge(task_table, data.frame(taxon = c('fia','bbs'), rv = c('alpha_richness','alpha_richness','beta_td_sorensen_pa','beta_td_sorensen_pa'), random_effect = rep(c('spatial','nonspatial'),each=4), predictor_group = 'null'), all = TRUE)

taxon <- task_table$taxon[task]
scale_km <- task_table$scale[task]
predictors <- task_table$predictor_group[task]
random_effect <- task_table$random_effect[task]
rv <- task_table$rv[task]

# Function used for model fitting.
source('/mnt/research/nasabio/code/fit_mv_mm.r')

options(mc.cores = 3)

# Load data for appropriate taxon
if (taxon == 'bbs') {
  load('/mnt/research/nasabio/temp/bbs_scaling_dat.RData')
  geodat <- bbsgeo
  biodat <- bbsbio
  siteid <- 'rteNo'

  # Logit transform beta td
  biodat$beta_td_sorensen_pa <- qlogis(biodat$beta_td_sorensen_pa)

} else {
  load('/mnt/research/nasabio/temp/fia_scaling_dat.RData')
  geodat <- fiageo
  biodat <- fiabio
  siteid <- 'PLT_CN'
  
  # Logit transform beta td.
  biodat$beta_td_sorensen_pa <- qlogis(biodat$beta_td_sorensen_pa)
  biodat$beta_td_sorensen <- qlogis(biodat$beta_td_sorensen)
  
  # Edit 18 June: Remove one unexplained outlier (~5E14) of elevation TRI
  geodat$elevation_30m_tri_100_mean[geodat$elevation_30m_tri_100_mean > 1e6] <- NA
}

# Get correct subset of predictors
prednames <- names(geodat)[!(names(geodat) %in% c('rteNo','HUC4','TNC','BCR'))]
prednames <- grep(scale_km, prednames, value = TRUE)

if (predictors == 'climate') prednames <- grep(paste(climate_prednames, collapse='|'), prednames, value = TRUE)
if (predictors == 'geo') prednames <- grep(paste(geo_prednames, collapse='|'), prednames, value = TRUE)
if (predictors == 'human') prednames <- grep(paste(human_prednames, collapse='|'), prednames, value = TRUE)
if (predictors == 'null') prednames <- character(0)

distrib <- 'gaussian' # Fit all with Gaussian
ecoregion <- 'TNC' # Use TNC as the ecoregion

added_priors <- NULL # For most tasks, model will converge using only default priors
# Edit 18 June: Add priors for the models that did not converge
if (rv == 'alpha_richness' & taxon == 'fia') {
  added_priors <- c(brms::set_prior('student_t(5, 0, 2)', class = 'Intercept'))
} 
if (rv == 'beta_td_sorensen_pa' & taxon == 'bbs') {
  added_priors <- c(brms::set_prior('lognormal(1, 1)', class = 'sdcar'))
}


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
				  random_effect_type = random_effect
				  )

# Save all fits
save(fit, file = paste0('/mnt/research/nasabio/temp/scalingfits/fit',task,'.RData'))