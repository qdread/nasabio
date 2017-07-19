# Fit occupancy model with same toy dataset in stan instead of jags.

rm(list = ls())
library(dplyr)
library(rstan)

load('C:/Users/Q/Dropbox/projects/nasabiodiv/occmod_workspace.r')
occ_model <- stan_model(file='occupancy/occmod.stan')

yr <- 1997 # Test year: 1997

###Combine input data for 1997: detection array and covariate information
X <- bbs_arrays$x[[which(bbs_arrays$year == yr)]]
X[is.na(X)] <- 0
elev <- bbselevvectors$x[[which(bbselevvectors$year == yr)]]

# Take a random sample of 50 sites to use for test data, to speed up the model fitting
sites_to_sample <- 50

set.seed(5156)
elev_subset <- sample(elev[!is.na(elev)], sites_to_sample)
X_subset <- X[names(elev_subset), ,]

sp.data <- list(nspec=dim(X_subset)[3], nsite=dim(X_subset)[1], nrep=dim(X_subset)[2], npred=1, X=X_subset, preds=matrix(as.numeric(scale(elev_subset)), ncol = 1))

# Here, no attempt is made to set reasonable initial values. If it doesn't work well, we can change that.

# Test sampling a very small number of iterations. 
occ_stanfit <- sampling(occ_model, data = sp.data, chains = 1, iter = 200, warmup = 100, thin = 1, init = 0.1)
# It didn't give any rejection messages, so it seems like the 

# Test sampling a small number of iterations with a very big dataset.

yr <- 2005

###Combine input data for 1997: detection array and covariate information
X <- bbs_arrays$x[[which(bbs_arrays$year == yr)]]
X[is.na(X)] <- 0
elev <- bbselevvectors$x[[which(bbselevvectors$year == yr)]]

sp.data2005 <- list(nspec=dim(X)[3], nsite=dim(X)[1], nrep=dim(X)[2], npred=1, X=X, preds=matrix(as.numeric(scale(elev)), ncol = 1))

occ_stanfit2005 <- sampling(occ_model, data = sp.data2005, chains = 1, iter = 200, warmup = 100, thin = 1, init = 0.1)

# Test sampling "random-effect only model" with a lot of data.
occ_ranef <- stan_model('occupancy/occmod_ranef.stan')
data_ranef2005 <- list(nspec=dim(X)[3], nsite=dim(X)[1], nrep=dim(X)[2], X=X)

occ_ranef_fit2005 <- sampling(occ_ranef, data = data_ranef2005, chains = 1, iter = 200, warmup = 100, thin = 1, init = 0.1)

# Dump the 2005 data to the hpcc
with(data_ranef2005, stan_rdump(names(data_ranef2005), 'X:/code/occupancy/data_ranef2005.txt'))

### cmdstan syntax to sample model.
./occmod_ranef sample num_samples=100 num_warmup=100 thin=1 data file=./data_ranef2005.txt init=0.1 output file=./test_ranef_output.csv
