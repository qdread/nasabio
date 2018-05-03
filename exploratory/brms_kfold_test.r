# Test brms kfold

fakedat <- data.frame(x = 1:50, y = 5 * 1:50 + 5 + rnorm(10,0,3))

library(brms)
fakemod <- brm(y ~ x, data = fakedat)

fakek <- kfold(fakemod, save_fits = TRUE)

names(fakek)

lapply(fakek$fits[,1], predict, newdata = fakedat) # This gives the predicted values of the 90% each time that are 

library(purrr)
oos_pred <- map2(fakek$fits[,'fit'], fakek$fits[,'omitted'], function(fit, idx) {
  pred <- predict(fit, newdata = fakedat[idx,])
  as.data.frame(cbind(idx = idx, fakedat[idx,], pred))
})

# rmse for each fold

oos_rmse <- map_dbl(oos_pred, function(x) {
  sqrt(mean((x$y - x$Estimate)^2))
})
sqrt(sum(oos_rmse^2)/length(oos_rmse))

pred_full <- predict(fakemod)
sqrt(mean((fakedat$y - pred_full[,'Estimate'])^2))

# General function to do k fold and get RMSE from the fit

kfold_rmse <- function(fit, k) {
  require(brms)
  require(purrr)
  
  kf <- kfold(fit, K = k, chains = 2, cores = 2, iter = 2000, warmup = 1000, save_fits = TRUE)
  
  # Predicted values for the 10% not in each fold.
  oos_pred <- map2(kf$fits[,'fit'], kf$fits[,'omitted'], function(fit, idx) {
    pred <- predict(fit, newdata = fit$data[idx,])
    as.data.frame(cbind(idx = idx, y = fit$data[idx, 1], pred))
  })
  
  # RMSE for each fold
  oos_rmse <- map_dbl(oos_pred, function(x) {
    sqrt(mean((x$y - x$Estimate)^2))
  })
  oos_rmse_total <- sqrt(sum(oos_rmse^2)/k)
  
  return(list(kfold_estimates = kf$estimates, rmse_total = oos_rmse_total, rmse_fold = oos_rmse))
  
}

# Count up how many sites are in each region

load('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_spatial_mm_dat.RData')

bbs_n_bcr <- bbsgeo %>% group_by(BCR) %>% summarize(n = n())
bbs_n_tnc <- bbsgeo %>% group_by(TNC) %>% summarize(n = n()) # Three have only one 

filter(bbs_n_tnc, n<5) # All but 3 have at least five.

bbs_n_huc <- bbsgeo %>% group_by(HUC4) %>% summarize(n = n()) # 17 have less than 5. We could get rid of them?

load('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_spatial_mm_dat.RData')
fia_n_bcr <- fiageo %>% group_by(BCR) %>% summarize(n = n())
fia_n_tnc <- fiageo %>% group_by(TNC) %>% summarize(n = n()) 
fia_n_huc <- fiageo %>% group_by(HUC4) %>% summarize(n = n()) 

sort(fia_n_tnc$n) # Two have less than 5
sort(fia_n_huc$n) # Five have less than 5.



