# Code to define and fit the phylogenetic/spatial imputation model in Stan.

# Step 1: compile stan model ----------------------------------------------


library(rstan)

trait_model <- stan_model(file = 'imputation/phylo_spatial_trait.stan') # Read from GitHub.
traitmissing_model <- stan_model(file = 'imputation/phylo_spatial_missing.stan') # Throws a warning but still works!

# Step 2: Create list containing m,n,N,p,Y,X,Z,and R for model fit --------

# Source script to load all the trait and predictor data from the nasabio space on the hpcc
source('imputation/loadtraitsandpredictors.r')

# Source functions to make Stan data from GitHub
source('imputation/makestandata.r')

trait_data_list <- make_standatalist(traitdat = species_traits, predictors = species_covariates, phy = species_phylogeny, evolution_model = 'ou')

# Missing values test: two missing
species_traits_missing <- species_traits
species_traits_missing[1, 5] <- NA
species_traits_missing[3, 3] <- NA

trait_data_list_test <- make_standatalist(traitdat = species_traits[1:5,], predictors = species_covariates[1:5,1:3], phy = drop.tip(species_phylogeny, species_phylogeny$tip.label[!species_phylogeny$tip.label %in% species_traits$species[1:5]]), evolution_model = 'ou')

trait_data_list_missingtest <- make_standatalist_missing(traitdat = species_traits_missing[1:5,], predictors = species_covariates[1:5,1:3], phy = drop.tip(species_phylogeny, species_phylogeny$tip.label[!species_phylogeny$tip.label %in% species_traits$species[1:5]]), evolution_model = 'ou')


# Step 3: Set model fitting options ---------------------------------------

# These numbers can be changed. If convergence is happening easily, we can reduce the number of chains and the number of iterations. 
# If it is not converging, we can increase the number of iterations.
# If that still doesn't help with convergence, we will have to change the prior specifications inside the model.
# If that STILL doesn't work, we will have to change the actual parameterization of the model (let's hope that isn't necessary).

n_cores <- parallel::detectCores() # 4 on Q's machine

rstan_options(auto_write = TRUE)
options(mc.cores = n_cores) 

### set to very low for testing purposes. 
n_chains <- 2
n_iter <- 1000 
n_warmup <- 100
n_thin <- 1

### Here is a better set of options when we fit the full model.
# n_chains <- 4
# n_iter <- 99999
# n_warmup <- 10000
# n_thin <- 10

# Here, we can write a function to initialize each chain with different values, because it's bad to initialize all chains with 1 for every parameter.

# Step 4: Fit model -------------------------------------------------------

trait_fit <- sampling(trait_model, data = trait_data_list_test, chains = n_chains, iter = n_iter, warmup = n_warmup, thin = n_thin, init = 1)

# Step 5: Examine diagnostic plots and metrics ----------------------------

# Extract summary information on parameters
summ_fit <- summary(trait_fit)

# Display the summary info for some of the parameters: beta only
summ_fit$summary[grep('beta', row.names(summ_fit$summary)),]

head(summ_fit$summary)


# Diagnostic plots to make sure the models converged
stan_diag(trait_fit)

# Must load bayesplot because ggplot2 is no longer compatible with rstan.
library(bayesplot)
draws <- as.array(trait_fit, pars="beta")
mcmc_trace(draws)



# Step 6: Fit model with some of Y missing --------------------------------

# This should impute the missing values in Y.

trait_fit_missing <- sampling(traitmissing_model, data = trait_data_list_missingtest, chains = 1, iter = 200, warmup = 100, thin = 1, init = 1) # works!

summ_missing <- summary(trait_fit_missing) # See the Y_mis parameters! These are the imputed missing values!
