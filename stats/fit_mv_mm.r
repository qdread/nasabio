# Function to fit multivariate mixed model
# Spatial model with ecoregion random effect
# QDR NASABioXGeo
# Originally written in April but this source code file created 14 June 2018

# Modified 02 May 2019: Exclude entire regions explicitly
# Modified 01 May 2019: add option to exclude some data points and estimate them with a missing data model -- this requires complete overhaul of code to create formula syntax. (and removal of rescor in formula)
# Modified 15 August 2018: add option to force intercept through zero
# Modified 18 June 2018: add drop = FALSE to scale so that it does not give weird output if only 1 response variable
# Modified 15 June 2018: correct scale() function to just return a numeric vector without attributes (also debug this a bit)

fit_mv_mm <- function(pred_df, resp_df, pred_vars, resp_vars, id_var, region_var, adj_matrix, distribution = 'gaussian', priors = NULL, n_chains = 2, n_iter = 2000, n_warmup = 1000, delta = 0.9, random_effect_type = 'spatial', force_zero_intercept = FALSE, missing_data = FALSE, exclude_locations = character(0)) {
  require(dplyr)
  require(brms)
  require(reshape2)
  # Build formula and data
  if (missing_data) {
	missing_df <- resp_df[,c(id_var, 'missing')]
  }
  id_df <- pred_df[, c(id_var, region_var)]
  resp_df <- resp_df[, c(id_var, resp_vars)]
  resp_df[,-1] <- scale(resp_df[,-1, drop = FALSE])
  resp_df[,-1] <- lapply(resp_df[,-1, drop = FALSE], as.numeric)
  names(id_df)[2] <- 'region' # Make sure the name of the region is called region so that the random effects will specify correctly.
  region_name <- 'region'
  
  if (missing_data) {
	resp_var_names <- paste(resp_vars, '| mi()')
  } else {
	resp_var_names <- resp_vars
  }
  # below, create full formula string for the model with predictors
  # create much simpler one if there aren't predictors (edited 13 June)
  if (length(pred_vars) > 0) {
	pred_df <- pred_df[, c(id_var, pred_vars)]
	pred_df[,-1] <- scale(pred_df[,-1, drop = FALSE])
	pred_df[,-1] <- lapply(pred_df[,-1, drop = FALSE], as.numeric)
	pred_var_names <- names(pred_df)[-1]
	fixed_effects <- paste(pred_var_names, collapse = '+')
	intercepts <- if (force_zero_intercept) '0' else paste('(1|', region_name, ')', sep = '')
	random_effects <- paste(c(intercepts, paste('(', pred_var_names, ' - 1|', region_name, ')', sep = '')), collapse = '+')
	formula_string <- mvbf(flist = paste(resp_var_names, '~', fixed_effects, '+', random_effects), rescor = FALSE)
	dat <- Reduce(left_join, list(id_df, resp_df, pred_df)) %>% filter(complete.cases(.))
  } else {
	intercepts <- if (force_zero_intercept) '0 +' else ''
	formula_string <- mvbf(flist = paste(resp_var_names, '~', intercepts, paste('(1|', region_name, ')', sep = '')), rescor = FALSE)
	dat <- Reduce(left_join, list(id_df, resp_df)) %>% filter(complete.cases(.))
  }
  
  # Added 1 May 2019: Change any missing data rows to NA
  if (missing_data) {
	dat <- dat %>% left_join(missing_df)
	dat[dat$missing, rv] <- NA
  }
  
  # Added 2 May 2018: get rid of any region that has less than 5 sites.
  dat <- dat %>% group_by(region) %>% filter(n() >= 5)
  
  # Added 2 May 2019: exclude entire regions explicitly if they are listed in the exclude_locations argument
  if (length(exclude_locations) > 0) {
	dat <- dat %>% filter(!grepl(paste(exclude_locations, collapse = '|'), region))
  }
  
  # Added 4 May 2018: if any region no longer has a neighbor at this point, get rid of it too.
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
