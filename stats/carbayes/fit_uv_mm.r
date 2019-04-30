# Function to fit univariate CAR model in CARBayes
# Spatial model with ecoregion random effect
# QDR NASABioXGeo

# Version created from fit_mv_mm.r on 30 April 2019
# Modified 15 August 2018: add option to force intercept through zero
# Modified 18 June 2018: add drop = FALSE to scale so that it does not give weird output if only 1 response variable
# Modified 15 June 2018: correct scale() function to just return a numeric vector without attributes (also debug this a bit)

fit_uv_mm <- function(pred_df, resp_df, pred_vars, resp_var, id_var, region_var, adj_matrix, distribution = 'gaussian', priors = NULL,n_iter = 2000, n_warmup = 1000, random_effect_type = 'spatial') {
  require(dplyr)
  require(CARBayes)
  require(reshape2)
  # Build formula and data
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_var)]
  resp_df[,-1] <- scale(resp_df[,-1, drop = FALSE])
  resp_df[,-1] <- lapply(resp_df[,-1, drop = FALSE], as.numeric)
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  # below, create full formula string for the model with predictors
  # create much simpler one if there aren't predictors (edited 13 June)
  if (length(pred_vars) > 0) {
	pred_df <- pred_df[, c(id_var, pred_vars)]
	pred_df[,-1] <- scale(pred_df[,-1, drop = FALSE])
	pred_df[,-1] <- lapply(pred_df[,-1, drop = FALSE], as.numeric)
	pred_var_names <- names(pred_df)[-1]
	fixed_effects <- paste(pred_var_names, collapse = '+')
	formula_string <- paste(resp_var, '~', fixed_effects)
	# Remove the NA covariate rows (not needed to do for response variable)
		
	dat <- Reduce(left_join, list(id_df, resp_df, pred_df))
  } else {
	formula_string <- paste(resp_var, '~', 1)
	dat <- Reduce(left_join, list(id_df, resp_df))
  }
    
  # Added 2 May: get rid of any region that has less than 5 sites.
  dat <- dat %>% group_by(region) %>% filter(n() >= 5)
  
  # Added 4 May: if any region no longer has a neighbor at this point, get rid of it too.
  reduced_adj_matrix <- adj_matrix[rownames(adj_matrix) %in% dat$region, rownames(adj_matrix) %in% dat$region]
  nneighb <- rowSums(reduced_adj_matrix)
  keep_regions <- names(nneighb)[nneighb > 0]
  dat <- filter(dat, region %in% keep_regions)
  
  # Fit model, extract coefficients, and format them
  if (random_effect_type == 'spatial') {
	rho <- NULL 
  } else {
	rho <- 0
  } 
  mm <- S.CARmultilevel(formula = formula_string, family = 'gaussian', data = dat, W = reduced_adj_matrix, ind.area = as.numeric(factor(dat$region)), burnin = NW, n.sample = NW + NI, thin = 1, rho = rho)
  
  # Extract fixed effects and random effects from the model, to save time later.
  
  # Edit 16 Aug: do not extract fixed effects (and combined fixed+random effects) if it is a null model without fixed effects.
  random_effects <- ranef(mm)
  random_effects <- cbind(effect = 'random', melt(random_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
  if (!force_zero_intercept | length(pred_vars) > 0) {
	fixed_effects <- fixef(mm)
    region_effects <- coef(mm)
	fixed_effects <- cbind(effect = 'fixed', region = as.character(NA), melt(fixed_effects, varnames = c('parameter', 'stat')))
	region_effects <- cbind(effect = 'coefficient', melt(region_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
    mm_coef <- fixed_effects %>% full_join(random_effects) %>% full_join(region_effects)
  } else {
	mm_coef <- random_effects
  }
  return(list(model = mm, coef = mm_coef))
}
