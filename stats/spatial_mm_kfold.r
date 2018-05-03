# Do K-fold cross validation on the spatial mixed models and return CV predictions and information criteria.
# QDR/NASABIOXGEO/02 May 2018

# Edited 3 May 2018: Make number of iterations and warmups variable. Also make the fold variable here.

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Pass in number of iterations as an argument.
args <- commandArgs(TRUE)

for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))
}

kfold_rmse <- function(fit, k, n_chains, n_iter, n_warmup, seed = 101) {
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
  
  kf <- kfold(fit, chains = n_chains, cores = n_chains, iter = n_iter, warmup = n_warmup, save_fits = TRUE, group = 'fold')
  
  # Predicted values for the subset not in each fold.
  oos_pred <- map2(kf$fits[,'fit'], kf$fits[,'omitted'], function(fit_fold, idx) {
    pred <- predict(fit_fold, newdata = fit$data[idx,])
    as.data.frame(cbind(idx = idx, y = fit$data[idx, 1], pred))
  })
  
  # RMSE for each fold
  oos_rmse <- map_dbl(oos_pred, function(x) {
    sqrt(mean((x$y - x$Estimate)^2))
  })
  oos_rmse_total <- sqrt(sum(oos_rmse^2)/k)
  
  # Put predicted values into a single data frame
  foldid <- rep(1:k, map_int(oos_pred, nrow))
  oos_pred <- data.frame(fold = foldid, do.call(rbind, oos_pred))
  oos_pred <- oos_pred[order(oos_pred$idx), ]
  
  return(list(kfold_estimates = kf$estimates, rmse_total = oos_rmse_total, rmse_fold = oos_rmse, oos_pred = oos_pred))
  
}

fp <- '/mnt/research/nasabio/temp/spammfit'

load(file.path(fp, paste0('fit', task, '.RData')))

kf <- kfold_rmse(fit$model, k = 5, n_chains = 2, n_iter = NI, n_warmup = NW, seed = task + 101)	
	
save(kf, file = paste0('/mnt/research/nasabio/temp/spammkfold/kfold_', task, '.RData'))
