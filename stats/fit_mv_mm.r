# Function to fit multivariate mixed model
# Spatial model with ecoregion random effect
# QDR NASABioXGeo
# Originally written in April but this source code file created 14 June 2018

# Modified 15 June 2018: correct scale() function to just return a numeric vector without attributes (also debug this a bit)

fit_mv_mm <- function(pred_df, resp_df, pred_vars, resp_vars, id_var, region_var, adj_matrix, distribution = 'gaussian', priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.9, random_effect_type = 'spatial') {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_vars)]
  resp_df[,-1] <- scale(resp_df[,-1])
  resp_df[,-1] <- lapply(resp_df[,-1], as.numeric)
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  resp_var_names <- paste0('cbind(', paste(resp_vars, collapse = ','), ')')
  # below, create full formula string for the model with predictors
  # create much simpler one if there aren't predictors (edited 13 June)
  if (length(pred_vars) > 0) {
	pred_df <- pred_df[, c(id_var, pred_vars)]
	pred_df[,-1] <- scale(pred_df[,-1])
	pred_df[,-1] <- lapply(pred_df[,-1], as.numeric)
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
  if (random_effect_type == 'spatial') {
	  mm <- brm(formula = formula_string, data = dat, family = distribution, autocor = cor_car(adj_matrix, formula = ~ 1|region),
				chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  } else {
	  mm <- brm(formula = formula_string, data = dat, family = distribution,
				chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  }
  fixed_effects <- fixef(mm)
  random_effects <- ranef(mm)
  region_effects <- coef(mm)
  fixed_effects <- cbind(effect = 'fixed', region = as.character(NA), melt(fixed_effects, varnames = c('parameter', 'stat')))
  random_effects <- cbind(effect = 'random', melt(random_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
  region_effects <- cbind(effect = 'coefficient', melt(region_effects$region, varnames = c('region', 'stat', 'parameter'))) %>% mutate(region = as.character(region))
  mm_coef <- fixed_effects %>% full_join(random_effects) %>% full_join(region_effects)
  return(list(model = mm, coef = mm_coef))
}
