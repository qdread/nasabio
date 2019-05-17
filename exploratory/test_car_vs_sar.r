# Simulation study demonstrating the bias-variance tradeoff with CAR and SAR models.
# QDR / nasabioxgeo / 10 May 2019

library(brms)
library(tidyverse)
options(mc.cores = 3)

# Load data and select a very small subset of it.
load('/mnt/research/nasabio/temp/bbs_spatial_mm_dat_50k.RData')

# Use only the ones touching Southern Blue Ridge and Piedmont
n_use <- grep('Blue Ridge|Piedmont', dimnames(tnc_bin)[[1]])[-3]
all_use <- apply(tnc_bin[n_use, ,drop=FALSE], 1, function(x) which(x==1))
all_use_sort <- sort(Reduce(union, all_use))

W <- tnc_bin[all_use_sort, all_use_sort]

# Use these 10 regions but only select 8 data points from each region

bbsdatasubset <- bbsgeo %>%
  filter(TNC %in% dimnames(W)[[1]]) %>%
  left_join(bbsbio, by = 'rteNo') %>%
  select(rteNo, TNC, alpha_richness, bio1_5k_50_mean, bio12_5k_50_mean) %>%
  group_by(TNC) %>%
  sample_n(8) %>%
  ungroup() %>%
  arrange(rteNo)

# Get the actual distances among the sites
bbscentroids <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv') %>%
  filter(rteNo %in% bbsdatasubset$rteNo) %>%
  arrange(rteNo)
bbs_dist <- dist(bbscentroids[,4:5]) %>% as.matrix
dimnames(bbs_dist)[[1]] <- bbscentroids$rteNo

fit_car <- brm(alpha_richness ~ bio1_5k_50_mean, data = bbsdatasubset, 
           family = 'gaussian', autocor = cor_car(W = W, formula = ~ 1|TNC), chains = 3) 
summary(fit_car)

fit_sar <- brm(alpha_richness ~ bio1_5k_50_mean, data = bbsdatasubset, 
               family = 'gaussian', autocor = cor_sar(bbs_dist, type = 'error'), chains = 3) 
summary(fit_sar)

# Do cross-validation for the models.
kfold_car <- kfold(fit_car, K = 5, chains = 1, save_fits = TRUE)
kfold_sar <- kfold(fit_sar, K = 5, chains = 1, save_fits = TRUE)

# Extract observed values, in-sample predicted values, and out-of-sample predicted values
insample_predict_car <- predict(fit_car, summary = FALSE)
insample_predict_sar <- predict(fit_sar, summary = FALSE)
kfold_predict_car <- kfold_predict(kfold_car)
kfold_predict_sar <- kfold_predict(kfold_sar)

# RMSEs
rmse <- function(y, yrep) {
  yrep_mean <- colMeans(yrep)
  sqrt(mean((yrep_mean - y)^2))
}


rmse_car <- with(kfold_predict_car, rmse(y, yrep))
rmse_sar <- with(kfold_predict_sar, rmse(y, yrep))


# Builtin example ---------------------------------------------------------

data(oldcol, package = "spdep")
fit1 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD, 
            autocor = cor_lagsar(COL.nb), 
            chains = 2, cores = 2)
summary(fit1)
plot(fit1)

library(spdep)

colnbmat <- nb2mat(COL.nb, style = 'B')
colnbmat[upper.tri(colnbmat)] <- t(colnbmat)[upper.tri(colnbmat)]

fit2 <- brm(CRIME ~ INC + HOVAL, data = COL.OLD, 
            autocor = cor_car(colnbmat), 
            chains = 2, cores = 2)

