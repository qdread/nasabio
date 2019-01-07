# R script to fit multivariate spatial mixed-effect models.
# New script forked from spatial_mm_parallel.r
# QDR/Nasabioxgeo/11 May 2018

# Modified 04 Jan 2019: edit file paths and task names for new OS
# Modified 25 July: tighten prior on intercepts for beta
# Modified 1 July: Geodiversity only as well as climate only
# Modified 14 June: Replace SD on the predictors with TRI
# Modified 14 June: Add priors to some of the beta-diversity models that didn't converge.
# Modified 13 June: Get rid of precipitation SD, include null models with space only and with space+climate only
# Modified 30 May: scale response variables in addition to predictors.
# Modified 29 May: add priors so that FIA alpha model can converge.
# Modified 14 May: take logit transformation of beta TD so that all can be modeled with multivariate normal.

# This time, only fit 50 km radius, get rid of human footprint since it does not fit with anything about geodiversity, and only use TNC.
# Also only use incidence-based for FIA.

task <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

NC <- 3
# If arguments are not specified, give default values
NI <- as.numeric(Sys.getenv('NI'))
NW <- as.numeric(Sys.getenv('NW'))
delta <- as.numeric(Sys.getenv('delta'))
if (is.na(NI)) NI <- 5000
if (is.na(NW)) NW <- 3000
if (is.na(delta)) delta <- 0.8

prednames <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')
climate_prednames <- c('bio1_5k_50_mean', 'bio12_5k_50_mean')
geo_prednames <- c('elevation_5k_tri_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'dhi_gpp_5k_tri_50_mean')
alpha_resp <- c('alpha_richness', 'alpha_phy_pa', 'alpha_func_pa')
beta_resp <- c('beta_td_sorensen_pa', 'beta_phy_pa', 'beta_func_pa')
gamma_resp <- c('gamma_richness', 'gamma_phy_pa', 'gamma_func_pa')

task_table <- data.frame(taxon = rep(c('fia','bbs'), each = 3),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
						 model = rep(c('full','climate','space', 'geo'),each=6),
                         stringsAsFactors = FALSE)

taxon <- task_table$taxon[task]
if(task_table$rv[task] == 'alpha') rv <- alpha_resp
if(task_table$rv[task] == 'beta') rv <- beta_resp
if(task_table$rv[task] == 'gamma') rv <- gamma_resp
if(task_table$model[task] == 'climate') prednames <- climate_prednames
if(task_table$model[task] == 'geo') prednames <- geo_prednames
if(task_table$model[task] == 'space') prednames <- character(0)

ecoregion <- task_table$ecoregion[task]

source('/mnt/research/nasabio/code/fit_mv_mm.r')

# Fit the model for the given response variable, taxon, and ecoregion
options(mc.cores = 3)

if (taxon == 'bbs') {
  load('/mnt/research/nasabio/temp/bbs_spatial_mm_dat_50k.RData')
  geodat <- bbsgeo
  biodat <- bbsbio
  siteid <- 'rteNo'

  # Added 14 May: logit transform beta td.
  biodat$beta_td_sorensen_pa <- qlogis(biodat$beta_td_sorensen_pa)

} else {
  load('/mnt/research/nasabio/temp/fia_spatial_mm_dat_50k.RData')
  geodat <- fiageo
  biodat <- fiabio
  siteid <- 'PLT_CN'
  # Added 14 May: logit transform beta td.
  biodat$beta_td_sorensen_pa <- qlogis(biodat$beta_td_sorensen_pa)
  biodat$beta_td_sorensen <- qlogis(biodat$beta_td_sorensen)
}


# Modified 14 May: model all with Gaussian
distrib <- 'gaussian'
         
# Priors (added May 29)
# --------------------

# Edit 04 Jan 2019: temporarily remove all priors (add some back in on 05 Jan)
# Edit May 31: Add priors for FIA intercepts and for BBS alpha sdcar
# Edit June 14: Add sdcar priors and intercept priors on FIA beta, sd car priors on BBS beta
library(brms, lib.loc = '/mnt/home/qdr/R/x86_64-pc-linux-gnu-library/3.5')
# 1st arg is df, 2nd is mu, 3rd is sigma for student t distribution
added_priors <- NULL
if (task_table$rv[task] == 'alpha' & taxon == 'fia') {
  added_priors <- c(set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alpharichness'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphaphypa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphafuncpa') )
} 
# if (task_table$rv[task] == 'beta' & taxon == 'fia') {
  # added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betatdsorensenpa'),
					# set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betaphypa'),
					# set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betafuncpa'),
					# set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betatdsorensenpa'),
					# set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betaphypa'),
					# set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betafuncpa') )					
# }
if (task_table$rv[task] == 'beta' & taxon == 'fia') {
  added_priors <- c(set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betatdsorensenpa'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betaphypa'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betafuncpa') )					
}
if (task_table$rv[task] == 'gamma' & taxon == 'fia') {
  added_priors <- c(set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammarichness'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammaphypa'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammafuncpa') )
} 
if (task_table$rv[task] == 'alpha' & taxon == 'bbs') {
  added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alpharichness'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alphaphypa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alphafuncpa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alpharichness'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphaphypa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphafuncpa') )
}
if (task_table$rv[task] == 'beta' & taxon == 'bbs') {
  added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betatdsorensenpa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betaphypa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betafuncpa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'betatdsorensenpa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'betaphypa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'betafuncpa') )
}
if (task_table$rv[task] == 'gamma' & taxon == 'bbs') {
  added_priors <- c(set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'gammarichness'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'gammaphypa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'gammafuncpa') )
} 
			 
# --------------------				  
				  
if (ecoregion == 'HUC4') eco_mat <- huc_bin
if (ecoregion == 'BCR') eco_mat <- bcr_bin
if (ecoregion == 'TNC') eco_mat <- tnc_bin

fit <- fit_mv_mm(pred_df = geodat, 
                      resp_df = biodat, 
                      pred_vars = prednames, 
                      resp_vars = rv, 
                      id_var = siteid, 
                      region_var = ecoregion, 
                      distribution = distrib, 
                      adj_matrix = eco_mat,
					  priors = added_priors,
                      n_chains = NC,
                      n_iter = NI,
                      n_warmup = NW,
					  delta = delta
					  )

# Save all fits
save(fit, file = paste0('/mnt/research/nasabio/temp/mvspam/fit',task,'.RData'))
