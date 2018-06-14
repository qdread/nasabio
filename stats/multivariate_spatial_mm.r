# R script to fit multivariate spatial mixed-effect models.
# New script forked from spatial_mm_parallel.r
# QDR/Nasabioxgeo/11 May 2018

# Modified 14 June: Add priors to some of the beta-diversity models that didn't converge.
# Modified 13 June: Get rid of precipitation SD, include null models with space only and with space+climate only
# Modified 30 May: scale response variables in addition to predictors.
# Modified 29 May: add priors so that FIA alpha model can converge.
# Modified 14 May: take logit transformation of beta TD so that all can be modeled with multivariate normal.

# This time, only fit 50 km radius, get rid of human footprint since it does not fit with anything about geodiversity, and only use TNC.
# Also only use incidence-based for FIA.

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

prednames <- c('elevation_5k_50_sd', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_50_sd')
climate_prednames <- c('bio1_5k_50_mean', 'bio12_5k_50_mean')
alpha_resp <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa")
beta_resp <- c("beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa")
gamma_resp <- c("gamma_richness", "gamma_phy_pa", "gamma_func_pa")

task_table <- data.frame(taxon = rep(c('fia','bbs'), each = 3),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
						 model = rep(c('full','climate','space'),each=6),
                         stringsAsFactors = FALSE)

taxon <- task_table$taxon[task]
if(task_table$rv[task] == 'alpha') rv <- alpha_resp
if(task_table$rv[task] == 'beta') rv <- beta_resp
if(task_table$rv[task] == 'gamma') rv <- gamma_resp
if(task_table$model[task] == 'climate') prednames <- climate_prednames
if(task_table$model[task] == 'space') prednames <- character(0)

ecoregion <- task_table$ecoregion[task]



fit_mv_mm <- function(pred_df, resp_df, pred_vars, resp_vars, id_var, region_var, adj_matrix, distribution = 'gaussian', priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.9) {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_vars)]
  resp_df[,-1] <- scale(resp_df[,-1])
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  resp_var_names <- paste0('cbind(', paste(resp_vars, collapse = ','), ')')
  # below, create full formula string for the model with predictors
  # create much simpler one if there aren't predictors (edited 13 June)
  if (length(pred_vars) > 0) {
	pred_df <- pred_df[, c(id_var, pred_vars)]
	pred_df[,-1] <- scale(pred_df[,-1])
	pred_var_names <- names(pred_df)[-1]
	fixed_effects <- paste(pred_var_names, collapse = '+')
	random_effects <- paste(c(paste('(1|', region_name, ')', sep = ''), paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
	formula_string <- paste(resp_var_names, '~', fixed_effects, '+', random_effects)
	dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  } else {
	formula_string <- paste(resp_var_names, '~', paste('(1|', region_name, ')', sep = ''))
	dat <- Reduce(left_join, list(id_df, resp_df)) %>% filter(complete.cases(.))
  }
    
  # Added 2 May: get rid of any region that has less than 5 sites.
  dat <- dat %>% group_by(region) %>% filter(n() >= 5)
  
  # Added 4 May: if any region no longer has a neighbor at this point, get rid of it too.
  reduced_adj_matrix <- adj_matrix[rownames(adj_matrix) %in% dat$region, rownames(adj_matrix) %in% dat$region]
  nneighb <- rowSums(reduced_adj_matrix)
  keep_regions <- names(nneighb)[nneighb > 0]
  dat <- filter(dat, region %in% keep_regions)
  
  # Fit model, extract coefficients, and format them
  mm <- brm(formula = formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region),
            chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  fixed_effects <- fixef(mm)
  random_effects <- ranef(mm)
  region_effects <- coef(mm)
  fixed_effects <- cbind(effect = 'fixed', region = as.character(NA), melt(fixed_effects, varnames = c('parameter', 'stat')))
  random_effects <- cbind(effect = 'random', melt(random_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
  region_effects <- cbind(effect = 'coefficient', melt(region_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
  mm_coef <- fixed_effects %>% full_join(random_effects) %>% full_join(region_effects)
  return(list(model = mm, coef = mm_coef))
}


# Fit the model for the given response variable, taxon, and ecoregion
options(mc.cores = 3)

if (taxon == 'bbs') {
  load('/mnt/research/nasabio/temp/bbs_spatial_mm_dat_50k.RData')
  geodat <- bbsgeo
  biodat <- bbsbio
  siteid <- 'rteNo'

  # Added April 30: Correction for outliers on beta functional
  # Edited June 4: don't get rid of those outliers.
  # biodat$beta_func_pa[biodat$beta_func_pa < -10] <- NA
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
# if (task_table$rv[task] == 'beta') {
	# distrib <- list('beta', 'gaussian', 'gaussian')
# } else {
	# distrib <- list('gaussian', 'gaussian', 'gaussian')
# }
                  
# Priors (added May 29)
# --------------------

# Edit May 31: Add priors for FIA intercepts and for BBS alpha sdcar
# Edit June 14: Add sdcar priors and intercept priors on FIA beta, sd car priors on BBS beta
library(brms)
# Tighten prior on the intercept for FIA alpha.
# 1st arg is df, 2nd is mu, 3rd is sigma for student t distribution
added_priors <- NULL
if (task_table$rv[task] == 'alpha' & taxon == 'fia') {
  added_priors <- c(set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alpharichness'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphaphypa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphafuncpa') )
} 
if (task_table$rv[task] == 'beta' & taxon == 'fia') {
  added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betatdsorensenpa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betaphypa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betafuncpa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'betatdsorensenpa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'betaphypa'),
					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'betafuncpa') )					
}
if (task_table$rv[task] == 'gamma' & taxon == 'fia') {
  added_priors <- c(set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammarichness'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammaphypa'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammafuncpa') )
} 
if (task_table$rv[task] == 'alpha' & taxon == 'bbs') {
  added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alpharichness'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alphaphypa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alphafuncpa') )
}
if (task_table$rv[task] == 'beta' & taxon == 'bbs') {
  added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betatdsorensenpa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betaphypa'),
					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betafuncpa') )
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
