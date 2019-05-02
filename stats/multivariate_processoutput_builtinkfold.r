# Summarize spatial multivariate mixed model output
# This version is for the models where I did the k-fold CV manually, built in to the main model script
# QDR/NASABIOXGEO/01 May 2019

# Modified 19 June: Get rid of the CV calculation (gives bad results)
# Modified 14 June: include newer "null" and subset models
# Modified 4 June: get RMSE for each iteration so we can put a credible interval on it.

task_table <- expand.grid(taxon = c('fia','bbs'),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
						 model = c('full','climate','space', 'geo'),
						 fold = 0:8,
                         stringsAsFactors = FALSE)

n_fits <- nrow(task_table)
n_full_fits <- sum(task_table$fold == 0)

fp <- '/mnt/research/nasabio/temp/mvspam'

library(brms)
library(purrr)
library(reshape2)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(abind)

registerDoParallel(cores = n_full_fits)

# ===============================================================================
# Inspect summaries. Only full fits, not the CV folds.
model_summ <- foreach (i = 1:n_full_fits) %dopar% {
  load(file.path(fp, paste0('fit', i, '.RData')))
  summary(fit$model)
}

# Check R-hat from summary.
model_rhat <- map(model_summ, function(x) {
  pars <- with(x, rbind(fixed, spec_pars, cor_pars, random[[1]]))
  pars[order(pars[,'Rhat'], decreasing = TRUE), ]
})

map(model_rhat, function(x) x[x[,'Rhat'] >= 1.1, ])
# ==============================================================================

# function to get correct response variable names and apply where needed.
get_correct_variable_names <- function(fit) {
	raw_resp_names <- fit$model$formula$responses
	resp_idx <- match(raw_resp_names, gsub('_', '', names(fit$model$data)))
	names(fit$model$data)[resp_idx]
}

# function to get relative RMSE for the model, including its MCMC quantiles
get_relative_rmse <- function(fit, resp_var_names, predicted_values) {
	# Prediction raw values. 
	pred_raw <- predict(fit$model, summary = FALSE)
	dimnames(pred_raw)[[3]] <- resp_var_names

	# Observed raw values
	obs_raw <- fit$model$data[, resp_var_names]

	# Get RMSE for each iteration and their quantiles
	rmse_quantiles <- sweep(pred_raw, 2:3, as.matrix(obs_raw), FUN = '-') %>% # Subtract predicted - observed
	melt(varnames = c('iter', 'idx', 'response')) %>%
	group_by(response, iter) %>%
	summarize(RMSE = sqrt(mean(value^2))) %>%
	ungroup %>% group_by(response) %>%
	summarize(RMSE_mean = sqrt(mean(RMSE^2)), 
			  RMSE_q025 = quantile(RMSE, probs = 0.025), 
			  RMSE_q975 = quantile(RMSE, probs = 0.975))

	# Generate ranges of observed data and divide this by the RMSE values to get the relative RMSE values
	predicted_values %>%
		group_by(response) %>%
		summarize(range_obs = diff(range(observed))) %>%
		left_join(rmse_quantiles) %>%
		mutate_at(vars(starts_with('RMSE')), funs(relative = ./range_obs))
}

model_stats <- foreach (i = 1:n_full_fits) %dopar% {
  load(file.path(fp, paste0('fit', i, '.RData')))
  
  # Get the correct variable names and apply them where needed.
  raw_resp_names <- fit$model$formula$responses
  resp_names <- get_correct_variable_names(fit)
  
  fit$coef <- fit$coef %>%
    separate(parameter, into = c('response', 'parameter'), extra = 'merge') %>%
    mutate(response = resp_names[match(response, raw_resp_names)])
  
  model_coef <- fit$coef # Includes fixed, random, and coefficient.
  model_pred <- predict(fit$model) # Returns raw array: n data x 4 stats x n response variables.
  # The predict() call takes a long time (~20 min or so in some cases)
  
  dimnames(model_pred)[[3]] <- resp_names
  model_pred <- melt(model_pred, varnames=c('idx','stat','response'))
  
  # Join predicted with observed values
  model_obs <- melt(cbind(idx = 1:nrow(fit$model$data), fit$model$data[, resp_names]), id.vars = 1, value.name = 'observed', variable.name = 'response')
  model_pred <- dcast(model_pred, idx + response ~ stat) %>%
    left_join(model_obs)
  
  # Here, do the RMSE for the model.
  model_rmse <- get_relative_rmse(fit, resp_names, model_pred)
  
  # Bayesian R-squared
  model_r2 <- cbind(task_table[i, ], response = resp_names, bayes_R2(fit$model))
  
    list(coef = model_coef, pred = model_pred, rmse = model_rmse, r2 = model_r2)

}

model_coef <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$coef)))
model_pred <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$pred)))
model_rmse <- map2_dfr(model_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], as.data.frame(x$rmse)))
model_r2 <- map_dfr(model_stats, 'r2')

write.csv(model_coef, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_coef.csv', row.names = FALSE)
write.csv(model_pred, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_pred.csv', row.names = FALSE)
write.csv(model_rmse, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_rmse.csv', row.names = FALSE)
write.csv(model_r2, '/mnt/research/nasabio/data/modelfits/multivariate_spatial_r2.csv', row.names = FALSE)

# K-fold output -----------------------------------------------------------

# Load data frame with the fold IDs in it, so we can see which to calculate the RMSE for.
fold_df <- read.csv('/mnt/research/nasabio/data/ecoregions/ecoregion_folds.csv', stringsAsFactors = FALSE)

# For each fit, load each k-fold subset (8) and combine the outputs.

# function to get the k-fold RMSE
get_kfold_rmse <- function(fit_ids, K = 8) {
	fit_folds <- map(fit_ids[-1], function(i) {
		load(file.path(fp, paste0('fit', i, '.RData')))
		fit
	})
	resp_names <- get_correct_variable_names(fit_folds[[1]])
	
	# Load full fit so we can get the data out
	load(file.path(fp, paste0('fit', fit_ids[1], '.RData')))
	
	pred_folds_raw <- map(fit_folds, function(x) {
		preds <- predict(x$model, summary = FALSE)
		dimnames(preds)[[3]] <- resp_names
		preds
	})
	
	# Extract slices of predicted and observed that correspond to the holdout data points for each fold.
	pred_folds_holdout <- map(1:K, function(i) {
		holdout_idx <- fit_folds[[i]]$model$data$region %in% fold_df$TNC[fold_df$fold == i]
		pred_folds_raw[[i]][, holdout_idx, ]
	})
	
	# also change this so it just reorders the main data df to the same order as the holdout
	obs_folds_holdout <- map(1:K, function(i) {
		holdout_idx <- fit$model$data$region %in% fold_df$TNC[fold_df$fold == i]
		fit$model$data[holdout_idx, c(resp_names)]
	})
	
	# Bind the slices into the proper dimensions 
	pred_all_holdout <- abind(pred_folds_holdout, along = 2)
	obs_all_holdout <- do.call('rbind', obs_folds_holdout)
	
	# sweep out observed from fitted and calculate RMSE
	rmse_quantiles <- sweep(pred_all_holdout, 2:3, as.matrix(obs_all_holdout), FUN = '-') %>% # Subtract predicted - observed
		melt(varnames = c('iter', 'idx', 'response')) %>%
		group_by(response, iter) %>%
		summarize(RMSE = sqrt(mean(value^2))) %>%
		ungroup %>% group_by(response) %>%
		summarize(RMSE_mean = sqrt(mean(RMSE^2)), 
				  RMSE_q025 = quantile(RMSE, probs = 0.025), 
				  RMSE_q975 = quantile(RMSE, probs = 0.975))

	# Generate ranges of observed data and divide this by the RMSE values to get the relative RMSE values
	obs_folds_holdout %>%
		melt(variable.name = 'response') %>%
		group_by(response) %>%
		summarize(range_obs = diff(range(value))) %>%
		left_join(rmse_quantiles) %>%
		mutate_at(vars(starts_with('RMSE')), funs(relative = ./range_obs))
	
}	

model_kfold_stats <- foreach(i = 1:n_full_fits) %dopar% {
	fit_ids <- with(task_table, which(taxon == taxon[i] & rv == rv[i] & model == model[i]))
	get_kfold_rmse(fit_ids, K = 8)
}

model_kfold_stats <- map2_dfr(model_kfold_stats, 1:n_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], x))
	
write.csv(model_kfold_stats, '/mnt/research/nasabio/data/modelfits/multivariate_kfold_rmse.csv', row.names = FALSE)
