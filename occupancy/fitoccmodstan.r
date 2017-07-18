# Fit occupancy model with same toy dataset in stan instead of jags.

rm(list = ls())
library(dplyr)
library(rstan)

load('C:/Users/Q/Dropbox/projects/nasabiodiv/occmod_workspace.r')
occ_model <- stan_model(file='occupancy/occmod.stan')

yr <- 1997 # Test year: 1997

###Combine input data for 1997: detection array and covariate information
X <- bbs_arrays$x[[which(bbs_arrays$year == 1997)]]
X[is.na(X)] <- 0
elev <- bbselevvectors$x[[which(bbselevvectors$year == 1997)]]

# Take a random sample of 50 sites to use for test data, to speed up the model fitting
sites_to_sample <- 50

set.seed(5156)
elev_subset <- sample(elev[!is.na(elev)], sites_to_sample)
X_subset <- X[names(elev_subset), ,]

sp.data <- list(nspec=dim(X_subset)[3], nsite=dim(X_subset)[1], nrep=dim(X_subset)[2], npred=1, X=X_subset, preds=matrix(as.numeric(scale(elev_subset)), ncol = 1))

# Here, no attempt is made to set reasonable initial values. If it doesn't work well, we can change that.

# Test sampling a very small number of iterations.
occ_stanfit <- sampling(occ_model, data = sp.data, chains = 1, iter = 200, warmup = 100, thin = 1, init = 0.1)
