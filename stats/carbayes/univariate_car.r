# R script to fit univariate CAR models with CARBayes
# QDR/Nasabioxgeo/30 April 2019

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

# If arguments are not specified, give default values
NI <- as.numeric(Sys.getenv('NI'))
NW <- as.numeric(Sys.getenv('NW'))
if (is.na(NI)) NI <- 5000
if (is.na(NW)) NW <- 3000

prednames <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')
climate_prednames <- c('bio1_5k_50_mean', 'bio12_5k_50_mean')
geo_prednames <- c('elevation_5k_tri_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'dhi_gpp_5k_tri_50_mean')
alpha_resp <- c('alpha_richness', 'alpha_phy_pa', 'alpha_func_pa')
beta_resp <- c('beta_td_sorensen_pa', 'beta_phy_pa', 'beta_func_pa')
gamma_resp <- c('gamma_richness', 'gamma_phy_pa', 'gamma_func_pa')

task_table <- expand.grid(taxon = c('fia', 'bbs'),
						  rv = c(alpha_resp, beta_resp, gamma_resp),
						  model = c('full','climate','space', 'geo'),
						  stringsAsFactors = FALSE)
						 
taxon <- task_table$taxon[task]
rv <- task_table$rv[task]
if(task_table$model[task] == 'climate') prednames <- climate_prednames
if(task_table$model[task] == 'geo') prednames <- geo_prednames
if(task_table$model[task] == 'space') prednames <- character(0)

source('/mnt/research/nasabio/code/carbayes/fit_uv_mm.r')

# Fit the model for the given response variable, taxon, and ecoregion

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

# For now don't use them
added_priors <- NULL
			 
# --------------------				  
				  
fit <- fit_uv_mm(pred_df = geodat, 
                      resp_df = biodat, 
                      pred_vars = prednames, 
                      resp_vars = rv, 
                      id_var = siteid, 
					  region_var = 'TNC',
                      distribution = distrib, 
                      adj_matrix = tnc_bin,
					  priors = added_priors,
                      n_iter = NI,
                      n_warmup = NW
					  )

# Save all fits
save(fit, file = paste0('/mnt/research/nasabio/temp/mvspam/uvfit',task,'.RData'))
