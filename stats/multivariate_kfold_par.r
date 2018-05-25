# Do K-fold cross validation on the multivariate mixed models and return CV predictions and information criteria.
# Using BRMS >2.3, this can be done with single subsets at a time in parallel.
# Parallel version created 25 May
# QDR/NASABIOXGEO/25 May 2018

# Total number of tasks is number of models * number of folds per model = 6 * 5 = 30
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
taskdf <- expand.grid(fold = 1:5, model = 1:6)
fit_n <- taskdf$model[task]
fold_n <- taskdf$fold[task]

# Pass in number of iterations as an argument.
args <- commandArgs(TRUE)

for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))
}

onefold <- function(fit, k, ksub, n_chains, n_iter, n_warmup, delta = 0.8, seed = 101) {
  require(brms)
  require(purrr)
  require(dplyr)
  
  # The group 'fold' is specified ahead of time.
  assign_fold <- function(n, k) {
    sample(rep_len(sample(1:k), n))
  }
  
  set.seed(seed)
  folds <- fit$data %>% group_by(region) %>% transmute(fold = assign_fold(n(), k))
  fit$data$fold <- factor(folds$fold)
  
  # Fit only the specified fold.
  kf <- kfold(fit, Ksub = as.array(ksub), chains = n_chains, cores = n_chains, iter = n_iter, warmup = n_warmup, control = list(adapt_delta = delta), save_fits = TRUE, group = 'fold')
  
  # Get names of response variables. Note that underscore characters are removed by brms, causing issues.
  resp_idx <- match(fit$formula$responses, gsub('_', '', names(fit$data)))
  resp_names <- names(fit$data)[resp_idx]
  
  # Predicted values for the subset not in the specified fold.
  # Since this is multivariate, we need to rewrite this code to get multiple y obs and y pred columns
  oos_pred <- map2(kf$fits[,'fit'], kf$fits[,'omitted'], function(fit_fold, idx) {
    pred <- predict(fit_fold, newdata = fit$data[idx,])
	# Convert each slice of pred into a list element (one for each response variable)
	lapply(1:(dim(pred)[3]), function(i) as.data.frame(cbind(idx = idx, y = fit$data[idx, resp_names[i]], pred[, , i])))
  })
  
  # RMSE for each response variable
  oos_rmse <- map_dbl(oos_pred$fit, function(x) {
    sqrt(mean((x$y - x$Estimate)^2))
  })
    
  # Add fold ID
  oos_pred <- data.frame(fold = ksub, yvar = rep(resp_names, map_int(oos_pred$fit, nrow)), do.call(rbind, oos_pred$fit))
  
  return(list(kfold_estimates = kf$estimates, rmse_fold = setNames(oos_rmse, resp_names), oos_pred = oos_pred))
  
}

fp <- '/mnt/research/nasabio/temp/mvspam'

load(file.path(fp, paste0('fit', fit_n, '.RData')))

kf <- onefold(fit$model, k = 5, ksub = fold_n, n_chains = 2, n_iter = NI, n_warmup = NW, delta = delta, seed = task + 303)	
	
save(kf, file = paste0('/mnt/research/nasabio/temp/mvspam/kfold_', fit_n, '_', fold_n, '.RData'))
