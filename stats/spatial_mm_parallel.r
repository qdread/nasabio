# R script to fit one spatial mixed model.
# Fit for a single combination of response variable x taxon (bbs/fia) x region (huc/bcr/tnc)
# QDR/Nasabioxgeo/25 Apr 2018

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

NC <- 3
NI <- 4000
NW <- 3000

prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')
respnames_fia <- c("alpha_richness", "alpha_effspn", "alpha_phy_pa",
                   "alpha_phy", "alpha_func_pa", "alpha_func", "beta_td_sorensen_pa",
                   "beta_td_sorensen", "beta_phy_pa", "beta_phy", "beta_func_pa",
                   "beta_func", "gamma_richness", "gamma_effspn", "gamma_phy_pa",
                   "gamma_phy", "gamma_func_pa", "gamma_func")
respnames_bbs <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa",
                   "beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa", "gamma_richness",
                   "gamma_phy_pa", "gamma_func_pa")

task_table <- data.frame(taxon = rep(c('fia','bbs'), c(length(respnames_fia), length(respnames_bbs))),
                         rv = c(respnames_fia, respnames_bbs),
                         ecoregion = rep(c('HUC4','BCR','TNC'), each = length(respnames_fia) + length(respnames_bbs)),
                         stringsAsFactors = FALSE)

taxon <- task_table$taxon[task]
rv <- task_table$rv[task]
ecoregion <- task_table$ecoregion[task]

fit_spatial_mm <- function(pred_df, resp_df, pred_vars, resp_var, id_var, region_var, adj_matrix, distribution = c('gaussian','beta'), n_chains = 2, n_iter = 2000, n_warmup = 1000) {
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


# Fit the model for the given response variable, taxon, and ecoregion
options(mc.cores = 3)

if (taxon == 'bbs') {
  load('/mnt/research/nasabio/temp/bbs_spatial_mm_dat.RData')
  geodat <- bbsgeo
  biodat <- bbsbio
  siteid <- 'rteNo'
} else {
  load('/mnt/research/nasabio/temp/fia_spatial_mm_dat.RData')
  geodat <- fiageo
  biodat <- fiabio
  siteid <- 'PLT_CN'
}

distrib <- ifelse(grepl('beta_td', rv), 'beta', 'gaussian') # Beta_td follows beta distribution, rest are normally distributed
                  
if (ecoregion == 'HUC4') eco_mat <- huc_bin
if (ecoregion == 'BCR') eco_mat <- bcr_bin
if (ecoregion == 'TNC') eco_mat <- tnc_bin

fit <- fit_spatial_mm(pred_df = geodat, 
                      resp_df = biodat, 
                      pred_vars = prednames, 
                      resp_var = rv, 
                      id_var = siteid, 
                      region_var = ecoregion, 
                      distribution = distrib, 
                      adj_matrix = eco_mat,
                      n_chains = NC,
                      n_iter = NI,
                      n_warmup = NW)

# Save all fits
save(fit, file = paste0('/mnt/research/nasabio/temp/spammfit/fit',task,'.RData'))
