# Run Jarzyna's occupancy model on our BBS data.
# Code modified by Q from the code Marta shared with me.
# Use the arrays and initial value matrices generated in a different script.
# Run separately for each year. Try to run in parallel if possible.
# Code to run rjags in parallel on the hpcc

# each year gets its own task

yr <- as.numeric(Sys.getenv('PBS_ARRAYID')) + 1996 # 20 tasks, 1997-2016

library(dplyr)
load('/mnt/research/nasabio/data/bbs/occmod_workspace.r')

###Combine input data for chosen year: detection array and covariate information
X <- bbs_arrays$x[[which(bbs_arrays$year == yr)]]
elev <- bbselevvectors$x[[which(bbselevvectors$year == yr)]]

# Create input data list, scaling the elevation predictor.
sp.data <- list(nspec=dim(X)[3], nsite=dim(X)[1], nrep=dim(X)[2], X=X, elev=as.numeric(scale(elev)))

###Specify the parameters to be monitored
sp.params <- c("Z", "occ_sp")
#Z matrix will store true prob of occurrence, so it is what we will be the most interested in
#additionally, occ_sp will store mean occupancy across all sites for all species

###Specify the initial values
# Latent variable Z, true occurrence, is initialized with the observed occurrences, which makes sense.
Zobs <- bbsbyroutelist$x[[which(bbsbyroutelist$year == yr)]]

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

# Set sampling options: number of burn-in steps, number of sampling steps, and how much to thin
# Use Jarzyna's settings
n_burn <- 500
n_thin <- 10
n_iter <- 20000
n_adapt <- 1000

library(runjags)
ocmod <- run.jags(model = modelString,
                  monitor = sp.params,
                  data = sp.data,
                  n.chains = 3,
                  inits = sp.inits(nspec = dim(X)[3], Zobs = Zobs),
                  burnin = n_burn,
                  sample = n_iter,
                  adapt = n_adapt,
                  thin = n_thin,
                  method = 'parallel')

# For now, just save the jags object
save(ocmod, file = paste0('/mnt/research/nasabio/data/bbs/jagsoutput/jagsout',yr,'.r'))