# Do K-fold cross validation on the spatial mixed models and return CV predictions and information criteria.
# QDR/NASABIOXGEO/02 May 2018

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

kfold_rmse <- function(fit, k) {
  require(brms)
  require(purrr)
  
  # The group 'fold' is specified ahead of time.
  kf <- kfold(fit, K = k, chains = 2, cores = 2, iter = 2000, warmup = 1000, save_fits = TRUE, group = 'fold')
  
  # Predicted values for the 10% not in each fold.
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
# For the ones that didn't converge the first time, load the model fit from a different directory
fp2 <- '/mnt/research/nasabio/temp/spammfit_moreiter'

# Run 10-fold cross validation on the fit.

if (!(file.exists(file.path(fp2, paste0('fit', task, '.RData'))))) {
  load(file.path(fp, paste0('fit', task, '.RData')))
} else {
  load(file.path(fp2, paste0('fit', task, '.RData')))
}
	
kf <- kfold_rmse(fit$model, k = 10)	
	
save(kf, file = paste0('/mnt/research/nasabio/temp/spammkfold/kfold_', task, '.RData'))
