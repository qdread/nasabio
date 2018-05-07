# Fit all the spatial models with no random effect.

NC <- 3
NI <- 5000
NW <- 3000
delta <- 0.8

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

prednames100 <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')
prednames50 <- c('elevation_5k_50_sd', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'bio12_5k_50_sd', 'dhi_gpp_5k_50_sd', 'human_footprint_5k_50_mean')
respnames_fia <- c("alpha_richness", "alpha_effspn", "alpha_phy_pa",
                   "alpha_phy", "alpha_func_pa", "alpha_func", "beta_td_sorensen_pa",
                   "beta_td_sorensen", "beta_phy_pa", "beta_phy", "beta_func_pa",
                   "beta_func", "gamma_richness", "gamma_effspn", "gamma_phy_pa",
                   "gamma_phy", "gamma_func_pa", "gamma_func")
respnames_bbs <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa",
                   "beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa", "gamma_richness",
                   "gamma_phy_pa", "gamma_func_pa")

# 27 models need to be fit for each radius. Probably not necessary to do parallel because it's pretty quick.
task_table <- data.frame(taxon = rep(c('fia','bbs'), c(length(respnames_fia), length(respnames_bbs))),
                         rv = c(respnames_fia, respnames_bbs),
						 radius = rep(c(50,100), each = length(respnames_fia)+length(respnames_bbs)),
                         stringsAsFactors = FALSE)

library(dplyr)
						 
task_table <- task_table %>%
	mutate(distrib = if_else(grepl('beta_td', rv), 'beta', 'gaussian'))

if (task == 1) task_table <- filter(task_table, radius == 50)
if (task == 2) task_table <- filter(task_table, radius == 100)	
	
fit_nospatial <- function(taxon, rv, radius, distrib, priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.8) {
  require(dplyr)
  require(brms)
  require(reshape2)
  
  # Get the correct data frames
  if (taxon == 'bbs' & radius == 100) {
	  load('/mnt/research/nasabio/temp/bbs_spatial_mm_dat.RData')
	  pred_df <- bbsgeo
	  resp_df <- bbsbio
	  id_var <- 'rteNo'
	  pred_vars <- 	prednames100
	  # Added April 30: Correction for outliers on beta functional
	  resp_df$beta_func_pa[resp_df$beta_func_pa < -10] <- NA

  } 
  
  if (taxon == 'fia' & radius == 100) {
	  load('/mnt/research/nasabio/temp/fia_spatial_mm_dat.RData')
	  pred_df <- fiageo
	  resp_df <- fiabio
	  id_var <- 'PLT_CN'
	  pred_vars <- 	prednames100
  }
  
  if (taxon == 'bbs' & radius == 50) {
	  load('/mnt/research/nasabio/temp/bbs_spatial_mm_dat_50k.RData')
	  pred_df <- bbsgeo
	  resp_df <- bbsbio
	  id_var <- 'rteNo'
	pred_vars <- 	prednames50
	  # Added April 30: Correction for outliers on beta functional
	  resp_df$beta_func_pa[resp_df$beta_func_pa < -10] <- NA

  } 
  
  if (taxon == 'fia' & radius == 50) {
	  load('/mnt/research/nasabio/temp/fia_spatial_mm_dat_50k.RData')
	  pred_df <- fiageo
	  resp_df <- fiabio
	  id_var <- 'PLT_CN'
	  pred_vars <- 	prednames50
  }
  
  # Build formula and data
  id_df <- pred_df[, id_var, drop = FALSE]
  resp_df <- resp_df[, c(id_var, rv)]
  pred_df <- pred_df[, c(id_var, pred_vars)]
  pred_df[,-1] <- scale(pred_df[,-1])
  pred_var_names <- names(pred_df)[-1]
  fixed_effects <- paste(pred_var_names, collapse = '+')
  formula_string <- paste(names(resp_df)[2], '~', fixed_effects)
  dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  
  # Fit model, extract coefficients, and format them
  mm <- brm(formula_string, data = dat, family = distrib, 
            chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  fixed_effects <- fixef(mm)
  fixed_effects <- cbind(effect = 'fixed', region = as.character(NA), melt(fixed_effects, varnames = c('parameter', 'stat')))
  mm_coef <- fixed_effects
  return(list(model = mm, coef = mm_coef))
}

library(purrr)
library(brms)
options(mc.cores=3)

fits <- pmap(task_table, fit_nospatial, 
									n_chains = NC,
									n_iter = NI,
									n_warmup = NW,
									delta = delta
			)
			
if (task == 1) file_name <- '/mnt/research/nasabio/temp/fits_nospatial_50k.RData'
if (task == 2) file_name <- '/mnt/research/nasabio/temp/fits_nospatial_100k.RData'
		
save(fits, file = file_name)

