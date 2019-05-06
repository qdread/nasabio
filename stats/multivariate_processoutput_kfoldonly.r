# Summarize spatial multivariate mixed model output
# Only K-fold, to be run remotely
# QDR/NASABIOXGEO/01 May 2019

# Modified 06 May 2019: extract predicted values within the first map() call so that fewer huge objects are created in workspace

K <- 63 # Number of regions with adequate number of data points

task_table <- expand.grid(taxon = c('fia','bbs'),
                         rv = c('alpha', 'beta', 'gamma'),
                         ecoregion = 'TNC',
						 model = c('full','climate','space', 'geo'),
						 fold = 0:K,
                         stringsAsFactors = FALSE)

n_fits <- nrow(task_table)
n_full_fits <- sum(task_table$fold == 0)

fp <- '/mnt/gs18/scratch/groups/nasabio/modelfits'

library(brms)
library(purrr)
library(reshape2)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(abind)

registerDoParallel(cores = 12)

# function to get correct response variable names and apply where needed.
get_correct_variable_names <- function(fit) {
	raw_resp_names <- fit$model$formula$responses
	resp_idx <- match(raw_resp_names, gsub('_', '', names(fit$model$data)))
	names(fit$model$data)[resp_idx]
}

# K-fold output -----------------------------------------------------------

# Load data frame with the fold IDs in it, so we can see which to calculate the RMSE for.
fold_df <- read.csv('/mnt/research/nasabio/data/ecoregions/ecoregion_folds.csv', stringsAsFactors = FALSE)
exclude_regions <- c('NA0801', 'NA0808', 'NA0417', 'NA0514', 'NA1202', 'NA1301')
region_folds <- fold_df$TNC
region_folds <- region_folds[!grepl(paste(exclude_regions, collapse = '|'), region_folds)]

# For each fit, load each k-fold subset (K = 63 for the regions) and combine the outputs.

# function to get the k-fold RMSE
get_kfold_rmse <- function(fit_ids, K) {
	pred_folds_holdout <- map(1:K, function(i) {
		# Load fit for fold i
		load(file.path(fp, paste0('fit', fit_ids[i+1], '.RData')))
		message('Fit ', fit_ids[i+1], ' loaded.')
		
		# Raw predicted values for fold i
		resp_names <- get_correct_variable_names(fit)
		preds <- predict(fit$model, summary = FALSE)
		dimnames(preds)[[3]] <- resp_names
		
		# Subset of fold i corresponding to holdout points
		holdout_idx <- fit$model$data$region %in% region_folds[i]
		
		message('Fit ', fit_ids[i+1], ' predicted.')
		
		preds[, holdout_idx, ]
	})
	
	# Load full fit so we can get the data out
	load(file.path(fp, paste0('fit', fit_ids[1], '.RData')))
	resp_names <- get_correct_variable_names(fit)
	
	# Get observed values corresponding to predicted values for calculating error.
	obs_folds_holdout <- map(1:K, function(i) {
		holdout_idx <- fit$model$data$region %in% region_folds[i]
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

model_kfold_stats <- foreach(i = 1:n_full_fits, .packages = c('brms', 'reshape2', 'reshape2', 'tidyr', 'dplyr', 'abind'), .export = c('get_correct_variable_names', 'region_folds')) %dopar% {
	fit_ids <- with(task_table, which(taxon == taxon[i] & rv == rv[i] & model == model[i]))
	kfold_rmse <- get_kfold_rmse(fit_ids, K = 63)
	message('Job ', i, ' done')
	kfold_rmse
}

model_kfold_stats <- map2_dfr(model_kfold_stats, 1:n_full_fits, function(x, y) cbind(taxon = task_table$taxon[y], rv = task_table$rv[y], ecoregion = task_table$ecoregion[y], model = task_table$model[y], x))
	
write.csv(model_kfold_stats, '/mnt/research/nasabio/data/modelfits/multivariate_kfold_rmse.csv', row.names = FALSE)
