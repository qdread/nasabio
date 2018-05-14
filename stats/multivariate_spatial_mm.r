# R script to fit multivariate spatial mixed-effect models.
# New script forked from spatial_mm_parallel.r
# QDR/Nasabioxgeo/11 May 2018

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

prednames <- c('elevation_5k_50_sd', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'bio12_5k_50_sd', 'dhi_gpp_5k_50_sd')
alpha_resp <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa")
beta_resp <- c("beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa")
gamma_resp <- c("gamma_richness", "gamma_phy_pa", "gamma_func_pa")

task_table <- data.frame(taxon = rep(c('fia','bbs'), each = 3),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
                         stringsAsFactors = FALSE)

taxon <- task_table$taxon[task]
if(task_table$rv[task] == 'alpha') rv <- alpha_resp
if(task_table$rv[task] == 'beta') rv <- beta_resp
if(task_table$rv[task] == 'gamma') rv <- gamma_resp
ecoregion <- task_table$ecoregion[task]



fit_mv_mm <- function(pred_df, resp_df, pred_vars, resp_vars, id_var, region_var, adj_matrix, distribution = 'gaussian', priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.9) {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_vars)]
  pred_df <- pred_df[, c(id_var, pred_vars)]
  pred_df[,-1] <- scale(pred_df[,-1])
  pred_var_names <- names(pred_df)[-1]
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  resp_var_names <- paste0('cbind(', paste(resp_vars, collapse = ','), ')')
  fixed_effects <- paste(pred_var_names, collapse = '+')
  random_effects <- paste(c(paste('(1|', region_name, ')', sep = ''), paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
  formula_string <- paste(resp_var_names, '~', fixed_effects, '+', random_effects)
  dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  
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
  biodat$beta_func_pa[biodat$beta_func_pa < -10] <- NA
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
                  
# Priors (added May 1)
# --------------------

# library(brms)
# added_priors <-  c(prior('normal(0,10)', class = 'b'),
				   # prior('lognormal(1,1)', class = 'sd', group = 'region'),
				   # prior('student_t(3,0,3)', class = 'sdcar')
				 # )
				 
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
					  #priors = added_priors,
                      n_chains = NC,
                      n_iter = NI,
                      n_warmup = NW,
					  delta = delta
					  )

# Save all fits
save(fit, file = paste0('/mnt/research/nasabio/temp/mvspam/fit',task,'.RData'))
