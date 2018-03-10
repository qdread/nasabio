# Fixed effects with confidence interval plots for mixed models
# Do in parallel because it is taking a long time
load('/mnt/research/nasabio/temp/mmfits.RData')

# Generate confidence intervals for fixed effects
library(lme4)
library(glmmTMB)
library(foreach)
library(doParallel)
registerDoParallel(cores = 9)

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

map_ci <- function(model_list) {
  foreach(i = model_list) %dopar% {
    model_ci(i)
  }
}

ci_bbs_bcr <- map_ci(fit_bbs_bcr)
ci_bbs_huc <- map_ci(fit_bbs_huc)
ci_bbs_tnc <- map_ci(fit_bbs_tnc)

ci_fia_bcr <- map_ci(fit_fia_bcr)
ci_fia_huc <- map_ci(fit_fia_huc)
ci_fia_tnc <- map_ci(fit_fia_tnc)

save(list = ls(pattern = 'ci_'), file = '/mnt/research/nasabio/temp/mmcis.RData')

# Create plots by region.
# X axis: predictor
# Y axis: coefficient for each predictor given diversity type (facet)
# Split FIA into incidence and abundance