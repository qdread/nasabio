# Do K-fold cross validation on the multivariate mixed models and return CV predictions and information criteria.
# Using BRMS >2.3, this can be done with single subsets at a time in parallel.
# Parallel version created 25 May
# QDR/NASABIOXGEO/25 May 2018

# New version created 29 Apr 2019: We are now doing SAR not CAR model because it is the only one that works with this type of CV.
# Modified 08 Apr 2019: change cross-validation scheme to use the preselected region chunks instead of doing stratified sampling. (10 chunks for 10-fold CV)
# Modified 05 Jan 2019: update for new OS.
# Modified 14 June: add more tasks. (again on 1 July)
# Modified 4 June: save raw predictions with each iteration so that we can put a credible interval on the k-fold RMSE too.

# Get arguments specified in sbatch
K <- as.numeric(Sys.getenv('K'))
NI <- as.numeric(Sys.getenv('NI'))
NW <- as.numeric(Sys.getenv('NW'))
delta <- as.numeric(Sys.getenv('delta'))

# Total number of tasks is number of models * number of folds per model = 24 * 10 = 240 ( or 24 * 8 = 192)
task <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
taskdf <- expand.grid(fold = 1:K, model = 1:24)
fit_n <- taskdf$model[task]
fold_n <- taskdf$fold[task]

library(brms, lib.loc = '/mnt/home/qdr/R/x86_64-pc-linux-gnu-library/3.5')

onefold <- function(fit, k, ksub, n_chains, n_iter, n_warmup, delta = 0.8, mcmc_seed = 101) {
  require(purrr)
  require(dplyr)
  require(reshape2)
  
  # Fit only the specified fold.
  kf <- kfold(fit, Ksub = as.array(ksub), chains = n_chains, cores = n_chains, iter = n_iter, warmup = n_warmup, control = list(adapt_delta = delta), save_fits = TRUE, group = 'fold', seed = mcmc_seed)
  
  # Get names of response variables. Note that underscore characters are removed by brms, causing issues.
  resp_idx <- match(fit$formula$responses, gsub('_', '', names(fit$data)))
  resp_names <- names(fit$data)[resp_idx]
  
  # Predicted values for the subset not in the specified fold.
  # Since this is multivariate, we need to rewrite this code to get multiple y obs and y pred columns
  oos_pred <- map2(kf$fits[,'fit'], kf$fits[,'omitted'], function(fit_fold, idx) {
    pred_raw <- predict(fit_fold, newdata = fit$data[idx,], summary = FALSE)
	dimnames(pred_raw)[[3]] <- resp_names
	obs_raw <- fit$data[idx, resp_names]
	sweep(pred_raw, 2:3, as.matrix(obs_raw), FUN = '-') %>% # Subtract predicted - observed
      melt(varnames = c('iter', 'idx', 'response'))
  })
  
  # Do not calculate the RMSE in here anymore. Do it once all the output has been combined so that we can get the credible interval.  

  # Add fold ID
  oos_pred <- data.frame(fold = ksub, bind_rows(oos_pred$fit))
  
  return(list(kfold_estimates = kf$estimates, oos_pred = oos_pred))
  
}

fp <- '/mnt/research/nasabio/temp/mvspam'

load(file.path(fp, paste0('fit', fit_n, '.RData')))

# Join fit data with ecoregion fold ID.
fold_df <- read.csv('/mnt/research/nasabio/data/ecoregions/ecoregion_folds.csv', stringsAsFactors = FALSE)
fold_byregion <- fold_df$fold[match(fit$model$data$region, fold_df$TNC)]
fit$model$data$fold <- factor(fold_byregion) 

kf <- onefold(fit$model, k = K, ksub = fold_n, n_chains = 2, n_iter = NI, n_warmup = NW, delta = delta, mcmc_seed = fit_n + fold_n + 303)	
	
save(kf, file = paste0('/mnt/research/nasabio/temp/mvspam/kfold_', fit_n, '_', fold_n, '.RData'))
