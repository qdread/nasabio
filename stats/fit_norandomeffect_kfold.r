# Do k-fold cross validation on fits (k=5)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

NI <- 5000
NW <- 3000
delta <- 0.8

if(task==1) load('/mnt/research/nasabio/temp/fits_nospatial_50k.RData')
if(task==2) load('/mnt/research/nasabio/temp/fits_nospatial_100k.RData')

kfold_rmse_nospatial <- function(fit, k, n_chains, n_iter, n_warmup, delta = 0.8, seed = 101) {
  require(brms)
  require(purrr)
  require(dplyr)
  
  set.seed(seed)
  
  kf <- kfold(fit, K = k, chains = n_chains, cores = n_chains, iter = n_iter, warmup = n_warmup, control = list(adapt_delta = delta), save_fits = TRUE)
  
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

library(purrr)
kfolds <- map2(fits, 10101 + 1:length(fits), function(x, y) kfold_rmse_nospatial(x$model, k = 5, n_chains = 2, n_iter = NI, n_warmup = NW, delta = delta, seed = y))

if(task==1) save(kfolds, file = '/mnt/research/nasabio/temp/kfolds_nospatial_50k.RData')
if(task==2) save(kfolds, file = '/mnt/research/nasabio/temp/kfolds_nospatial_100k.RData')
