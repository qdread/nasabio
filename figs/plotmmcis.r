# Fixed effects with confidence interval plots for mixed models
load('C:/Users/Q/Dropbox/projects/nasabiodiv/mmfits.RData')

# Generate confidence intervals for fixed effects
library(purrr)
library(lme4)
library(glmmTMB)

prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')

model_ci <- function(x,pnames=prednames) {
  if (inherits(x$model, 'lmerMod')) {
    ci <- confint(x$model, method = 'boot', nsim = 999)
    ci <- data.frame(pred = dimnames(ci)[[1]], q025 = ci[,1], q975 = ci[,2])
  } else {
    ci <- confint(x$model, method = 'wald')
    ci <- data.frame(pred = dimnames(ci)[[1]], q025 = ci[,1], q975 = ci[,2])
    ci$pred <- gsub('cond.', '', ci$pred)
  }
  ci[ci$pred %in% pnames, ]
}

ci_bbs_bcr <- map(fit_bbs_bcr, model_ci)
ci_bbs_huc <- map(fit_bbs_huc, model_ci)
ci_bbs_tnc <- map(fit_bbs_tnc, model_ci)


ci_fia_bcr <- map(fit_fia_bcr, model_ci)
ci_fia_huc <- map(fit_fia_huc, model_ci)
ci_fia_tnc <- map(fit_fia_tnc, model_ci)


# Create plots by region.
# X axis: predictor
# Y axis: coefficient for each predictor given diversity type (facet)
# Split FIA into incidence and abundance