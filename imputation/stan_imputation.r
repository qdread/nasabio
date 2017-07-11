# Code to define and fit the phylogenetic/spatial imputation model in Stan.

# Step 1: compile stan model ----------------------------------------------


library(rstan)

trait_model <- stan_model(file = 'imputation/phylo_spatial_trait.stan') # Read from GitHub.

# Step 2: Create list containing m,n,N,p,Y,X,Z,and R for model fit --------

# environmental predictors, phylogenetic VCV matrix, trait data, and spatial locations


# For now, create environmental predictors and traits from species means (mean locations)

# Use Stevens traits, and get median locations of each species from FIA dataset.
# Use the environmental predictors from FIA.

trait_training <- read.csv('X:/jay/tree_data/stevenstraits/trait_stevens_training.csv', stringsAsFactors = FALSE)

# Get median locations for each species from TRY
# All covariates for those locations have already been extracted from different GIS layers.
try_covariates <- read.csv('X:/data/fia/imputation_covariates.csv', stringsAsFactors = FALSE)

# Get only a subset of the covariates: elevation, minimum temperature, temperature seasonality, minimum precipitation, precipitation seasonality, lai, npp
try_covariates_reduced <- try_covariates %>%
  select(AccSpeciesName, lat_aea, lon_aea, elevation, bio6, bio4, bio14, bio15, lai, npp) %>%
  rename(temp_min_monthly = bio6, temp_seasonality = bio4, precip_min_monthly = bio14, precip_seasonality = bio15) %>%
  filter(complete.cases(.))


library(ape)
allfia_phylogeny <- read.tree('X:/data/fia/tree_all_final_031716.txt')

# Pare down the trait list and phylogeny to species that we have valid locations for. (64 spp)

species_list <- Reduce('intersect', list(try_covariates_reduced$AccSpeciesName, allfia_phylogeny$tip.label, trait_training$Scientific_Name))

species_covariates <- subset(try_covariates_reduced, AccSpeciesName %in% species_list)
species_phylogeny <- drop.tip(allfia_phylogeny, tip = allfia_phylogeny$tip.label[!allfia_phylogeny$tip.label %in% species_list])
species_traits <- subset(trait_training, Scientific_Name %in% species_list)

# Sort them all so that they are in the same order
species_covariates <- species_covariates %>%
  rename(species = AccSpeciesName) %>%
  arrange(species)

species_traits <- species_traits %>%
  rename(species = Scientific_Name) %>%
  arrange(species)

# Z is design matrix with dimensions (n trait * n indiv) x (n trait * n species) to make sure each individual maps to the right species

# Define function to make Z

make_Z <- function(traitdat) {
  Z <- matrix(0, nrow = nrow(traitdat), ncol = length(unique(traitdat$trait)) * length(unique(traitdat$species)))
  
  # For each row in traitdat, put a 1 in the correct row and column of Z.
  species_trait_combo <- with(traitdat, paste(species,trait))
  col_idx <- match(species_trait_combo, unique(species_trait_combo))
  for (i in 1:length(col_idx)) Z[i, col_idx[i]] <- 1
  return(Z)
}

# X is predictor matrix with dimensions (n trait * n indiv) x (n trait * n predictors)

# Define function to make X (Jay is writing)

make_X <- function(preds, m) {
  preds <- as.matrix(preds)
  p = ncol(preds) 
  n = nrow(preds)
  #X <- matrix(nrow=(m*nrow(preds)),ncol = m*(p+1), byrow = TRUE)
  bigList <- list()
  for(k in 1:n){
    mlist <- list()
    for(i in 1:m){
      
      X <- matrix(0,nrow=m,ncol=p+1) 
      X[i,] <- c(1,as.numeric(preds[k,]))
      mlist[[length(mlist) + 1]] <- X
    }
    ilist <- do.call(cbind,mlist)
    bigList[[k]] <- ilist
  }
  return(do.call(rbind,bigList))
}

# Y is vector of trait values by individuals/species, all strung into one long vector

# Define function to make Y
  
make_Y <- function(traitdat) {
  # traitdata has first column as species, rest as traits
  Y <- do.call(c, traitdat[,-1])
  names(Y) <- paste(traitdat$species, rep(names(traitdat)[-1], each = nrow(traitdat)), sep = '_')
  return(Y)
}
  
# Define function to make data list for stan model

make_standatalist <- function(traitdat, predictors, phy, evolution_model = c('brownian','ou')) {
  require(reshape2)
  require(Matrix)
  traitdat_melted <- melt(traitdat, id.vars = 'species', variable.name = 'trait')
  
  splist <- unique(traitdat$species)
  
  # R: phylogenetic VCV matrix. Can be either Brownian or Ornstein.
  if (evolution_model[1] == 'brownian') {
    R <- as.matrix(vcv(phy))
  }
  if (evolution_model[1] == 'ou') {
    # Function for vcv under OU model is in PIGShift
    library(PIGShift)
    R <- OU.vcv(phy, theta = 1) ### Need to check whether this is sensitive to our choice of theta.
    dimnames(R) <- list(phy$tip.label, phy$tip.label) # Label matrix rows and columns with species names.
  }
  
  R <- R[splist, splist] # Sort matrix so that it's in the same order as the trait and environment data frames.
  
  # Force R to be symmetric
  R <- as.matrix(forceSymmetric(R, "U"))
  
  m <- length(unique(traitdat_melted$trait))
  X <- make_X(predictors[,-1], m)
  Z <- make_Z(traitdat_melted)
  Y <- make_Y(traitdat)
  
  datlist <- list(m = m,
                  n = length(unique(traitdat$species)),
                  N = nrow(traitdat_melted)/m,
                  p = ncol(predictors), # this only works out because predictors has 1 extra column for species, then we add 1 to p for the intercept
                  X = X,
                  Y = Y,
                  Z = Z,
                  R = R)
  
  return(datlist)
}

trait_data_list <- make_standatalist(traitdat = species_traits, predictors = species_covariates, phy = species_phylogeny, evolution_model = 'ou')

trait_data_list_test <- make_standatalist(traitdat = species_traits[1:5,], predictors = species_covariates[1:5,1:3], phy = drop.tip(species_phylogeny, species_phylogeny$tip.label[!species_phylogeny$tip.label %in% species_traits$species[1:5]]), evolution_model = 'ou')

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
missing_data_test <- trait_data_list_test
missing_data_test$Y[1] <- NA

trait_fit_missing <- sampling(trait_model, data = missing_data_test, chains = 1, iter = 100, warmup = 100, thin = 1, init = 1) # FAILS
