# Function to fit multivariate mixed model
# SAR model
# QDR NASABioXGeo

# New version created 29 April 2019: This is now a SAR model not a CAR (continuous random effects not discrete based on regions) - note SAR does not support multivariate so we have to keep it to a single response variable.
# Modified 15 August 2018: add option to force intercept through zero
# Modified 18 June 2018: add drop = FALSE to scale so that it does not give weird output if only 1 response variable
# Modified 15 June 2018: correct scale() function to just return a numeric vector without attributes (also debug this a bit)

fit_mv_sar <- function(pred_df, resp_df, pred_vars, resp_vars, id_var, dist_matrix, distribution = 'gaussian', priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.9, random_effect_type = 'spatial', force_zero_intercept = FALSE) {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  resp_df <- resp_df[, c(id_var, resp_vars)]
  resp_df[,-1] <- scale(resp_df[,-1, drop = FALSE])
  resp_df[,-1] <- lapply(resp_df[,-1, drop = FALSE], as.numeric)
  resp_var_names <- resp_vars
  # below, create full formula string for the model with predictors
  # create much simpler one if there aren't predictors (edited 13 June)
  if (length(pred_vars) > 0) {
	pred_df <- pred_df[, c(id_var, pred_vars)]
	pred_df[,-1] <- scale(pred_df[,-1, drop = FALSE])
	pred_df[,-1] <- lapply(pred_df[,-1, drop = FALSE], as.numeric)
	pred_var_names <- names(pred_df)[-1]
	fixed_effects <- paste(pred_var_names, collapse = '+')
	#intercepts <- if (force_zero_intercept) '0' else paste('(1|', region_name, ')', sep = '')
	#random_effects <- paste(c(intercepts, paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
	#formula_string <- paste(resp_var_names, '~', fixed_effects, '+', random_effects)
	formula_string <- paste(resp_var_names, '~', fixed_effects)
	
  } else {
	#intercepts <- if (force_zero_intercept) '0 +' else ''
	#formula_string <- paste(resp_var_names, '~', intercepts, paste('(1|', region_name, ')', sep = ''))
	formula_string <- paste(resp_var_names, '~', 1)
  }
  
  dat <- left_join(resp_df, pred_df)
  use_obs <- complete.cases(dat)
  dat <- dat[use_obs,]
  dist_matrix <- dist_matrix[use_obs, use_obs]
    
  # Fit model, extract coefficients, and format them
  if (random_effect_type == 'spatial') {
	  mm <- brm(formula = formula_string, data = dat, family = distribution, autocor = cor_sar(dist_matrix, 'lag'),
				chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  } else {
	  mm <- brm(formula = formula_string, data = dat, family = distribution,
				chains = n_chains, iter = n_iter, warmup = n_warmup, prior = priors, control = list(adapt_delta = delta))
  }
  # Edit 16 Aug: do not extract fixed effects (and combined fixed+random effects) if it is a null model without fixed effects.
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
