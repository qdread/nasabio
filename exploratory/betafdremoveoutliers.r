# Try to remove outliers and see if beta FD model fit is any better


prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')

# Fit all BBS response variables to same predictors by HUC, BCR, TNC

library(purrr)

distribs <- ifelse(grepl('beta_td', names(bbsbio)[-1]), 'beta', 'normal') # Beta_td follows beta distribution, rest are normally distributed

fit_bbs_huc <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'HUC4', distr))
fit_bbs_bcr <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'BCR', distr))
fit_bbs_tnc <- map2(names(bbsbio)[-1], distribs, function(varname, distr) fit_mm(bbsgeo, bbsbio, prednames, varname, 'rteNo', 'TNC', distr))

fit_bcr_all <- fit_mm(bbsgeo, bbsbio, prednames, names(bbsbio)[9], 'rteNo', 'BCR', 'normal')

# Remove outliers
bbsbio_no_outlier <-  bbsbio %>% mutate(beta_fd_pairwise_pa_z_100 = if_else(beta_fd_pairwise_pa_z_100 < -25, as.numeric(NA), beta_fd_pairwise_pa_z_100))
fit_bcr_50 <- fit_mm(bbsgeo, bbsbio_no_outlier, prednames, names(bbsbio)[9], 'rteNo', 'BCR', 'normal')

