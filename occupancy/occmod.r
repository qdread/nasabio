# Run Jarzyna's occupancy model on our BBS data.
# Code modified by Q from the code Marta shared with me.
# Use the arrays and initial value matrices generated in a different script.
# Run separately for each year. Try to run in parallel if possible.

rm(list = ls())
library(dplyr)
load('C:/Users/Q/Dropbox/projects/nasabiodiv/occmod_workspace.r')

yr <- 1997 # Test year: 1997

###Combine input data for 1997: detection array and covariate information
X <- bbs_arrays$x[[which(bbs_arrays$year == 1997)]]
elev <- bbselevvectors$x[[which(bbselevvectors$year == 1997)]]

# Take a random sample of 50 sites to use for test data, to speed up the model fitting
sites_to_sample <- 50

set.seed(5156)
elev_subset <- sample(elev[!is.na(elev)], sites_to_sample)
X_subset <- X[names(elev_subset), ,]

# Create input data list, scaling the elevation predictor.
sp.data <- list(nspec=dim(X_subset)[3], nsite=dim(X_subset)[1], nrep=dim(X_subset)[2], X=X_subset, elev=as.numeric(scale(elev_subset)))

###Specify the parameters to be monitored
sp.params <- c("Z", "occ_sp")
#Z matrix will store true prob of occurrence, so it is what we will be the most interested in
#additionally, occ_sp will store mean occupancy across all sites for all species

###Specify the initial values
# Latent variable Z, true occurrence, is initialized with the observed occurrences, which makes sense.
Zobs <- bbsbyroutelist$x[[which(bbsbyroutelist$year == 1997)]]
Zobs_subset <- Zobs[names(elev_subset), ]

# Function to specify initial values.
sp.inits <- function(nspec, Zobs) {
  list(psi.mean=runif(1,0.001,0.99), 
       theta.mean=runif(1,0.001,0.99),
       u=rnorm(nspec), 
       v=rnorm(nspec),
       Z = Zobs,
       alpha1=rnorm(nspec))
}

#JAGS model
modelString <- "
    model {
    
    #Prior distributions on the community level occupancy and detection covariates
    psi.mean ~ dunif(0.001,0.99) #vague prior for the hyperparameter of the community-level occupancy covariates
    a <- log(psi.mean) - log(1-psi.mean)
    
    theta.mean ~ dunif(0.001,0.99) #vague prior for the hyperparameter of the community-level detection covariates
    b <- log(theta.mean) - log(1-theta.mean)
    mu.alpha1 ~ dnorm(0, 0.01) #vague Normal prior for of the community-level habitat (alpahs) and sampling (betas) covariates
    
    
    tau1 ~ dgamma(10,1)
    tau2 ~ dgamma(10,1)
    
    tau.alpha1 ~ dgamma(10,1)
    
    rho ~ dunif(-0.99,0.99)
    var.v <- tau2 /(1.-pow(rho,2))
    
    sigma1 <- 1/sqrt(tau1)
    sigma2 <- 1/sqrt(tau2)
    
    for (i in 1:nspec) {
    
    #Prior distributions for the occupancy and detection covariates for each species 
    u[i] ~ dnorm(a, tau1)
    
    mu.v[i] <- b + (rho*sigma2 /sigma1)*(u[i]-a)
    v[i] ~ dnorm(mu.v[i], var.v)
    
    alpha1[i] ~ dnorm(mu.alpha1, tau.alpha1)
    
    
    #Estimate the occupancy probability (latent Z matrix) for each species 
    #at each point (i.e., route or site)
    for (j in 1:nsite) {
    logit(psi[j,i]) <- u[i] + alpha1[i]*elev[j] 
    mu.psi[j,i] <- psi[j,i]
    Z[j,i] ~ dbin(psi[j,i], 1)#Z is generally not observed with certainty, instead
    #we observed data theta[i,j,k] for species i at site j during sampling period k
    
    #Estimate the species specific detection probability for every rep at each point where the species occurs (Z=1)
    # Q edited line 85 to say nrep instead of nrep[j]
    for (k in 1:nrep) {  
    logit(theta[j,k,i]) <- v[i] 
    mu.theta[j,k,i] <- theta[j,k,i]*Z[j,i]
    X[j,k,i] ~ dbin(mu.theta[j,k,i], 1) #X is the 3D array of dependent variable: The detection/non-
    #detection data is defined in a three dimensional
    #array X where the first dimension, j, is the point; the second
    #dimension, k, is the rep; and the last dimension, i, is the species.
    
    } } }
    
    #Estimation of species occupancy (averaged across all the sites)
    for(i in 1:nspec) {
    occ_sp[i] <- sum(Z[1:nsite,i])/nsite
    }
    
    #End model specification
    }
    "

# Test the model without parallel computing.
library(rjags)

# Set sampling options: number of burn-in steps, number of sampling steps, and how much to thin
# Use Jarzyna's settings
n_burn <- 500
n_thin <- 10
n_iter <- 20000

# Compile and initialize model
ocmod <- jags.model(file = textConnection(modelString), 
                    inits = sp.inits(nspec = dim(X_subset)[3], Zobs = Zobs_subset), 
                    data = sp.data, n.chains = 3)

# Do MCMC sampling including burn-in steps
# Even for the toy dataset with 10 sites, this takes a long time (non-parallel, roughly 30 min?)
# For 50 sites, it takes 2 or 3 hours.
update(ocmod, n.iter = n_burn)
out <- coda.samples(ocmod, n.iter = n_iter, variable.names = sp.params, thin = n_thin)

# Burn-in already removed, so concatenate each of the 3 chains to get the posterior.
out.mcmc <- as.mcmc(do.call(rbind, out))
str(out.mcmc)

# Create data frame with summary statistics for the output.
# There are 6000 posterior samples of 30753 parameters here
# 30753 = 603 species * (50 sites + 1)
# each species by site has an occupancy estimate (Z) plus a mean occupancy across sites (occ_sp), hence the 50+1.

summ_stats <- function(x) {
  res <- c(mean(x), sd(x), quantile(x, probs = c(0.025, 0.5, 0.975)))
  names(res) <- c('mean','sd','q025','median','q975')
  return(res)
}
 
output <- t(apply(out.mcmc, 2, summ_stats))

save(output, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/jagsoutput_occ_toydata.r')
save(Zobs_subset, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/jagsoutput_observations_toydata.r')
# Compare the observed values with the occupancy values

# Function to convert the MCMC summary output into a matrix with the same dimensions as the observed data
mcout2matrix <- function(mc) {
  library(stringr)
  dn <- dimnames(mc)[[1]]
  varn <- unlist(str_extract_all(dn, "[_[:alpha:]]+")) # letters or the underscore character.
  nums <- str_extract_all(dn, "[0-9]+")
  rowidx <- sapply(nums, function(x) as.integer(x[1]))
  colidx <- sapply(nums, function(x) as.integer(x[2]))
  n_sites <- max(rowidx)
  n_spp <- max(colidx, na.rm=TRUE)
  Z <- matrix(0, nrow = n_sites, ncol = n_spp)
  occ_sp <- numeric(n_spp)
  for (i in 1:nrow(mc)) {
    if (varn[i] == 'Z') Z[rowidx[i], colidx[i]] <- output[i, 'median']
    if (varn[i] == 'occ_sp') occ_sp[rowidx[i]] <- output[i, 'median']
  }
  return(list(Z=Z, occ_sp=occ_sp))
}

outmat <- mcout2matrix(output)

# Compare mean occupancy across the ten sites with mean observed presence.
meanobs <- apply(Zobs_subset, 2, function(x) sum(x)/length(x))

obspreddat <- data.frame(obs = meanobs, pred = outmat$occ_sp)

table(obspreddat$obs == obspreddat$pred)

library(cowplot)
ggplot(obspreddat, aes(x=obs, y=pred, color=obs==pred)) +
  geom_jitter() +
  scale_color_discrete(name = 'Any change?', labels=c('no (507 species)','yes (96 species)')) +
  theme(legend.position=c(0.2,0.8))

table(obspreddat$obs, obspreddat$pred)

# Histogram of changes in occupancy.
ggplot(obspreddat, aes(x = pred - obs)) + geom_histogram(binwidth = 0.02) +
  scale_y_continuous(expand = c(0,0))
