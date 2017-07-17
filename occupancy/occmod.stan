// Jarzyna's occupancy model translated from JAGS to STAN
// QDR 17 July 2017
// Goal is to speed up the fitting process (vectorized instead of with loops)
// Also generalize so that we can include multiple predictors.

data {
	int<lower=0> nspec;									// Number of species
	int<lower=0> nsite;									// Number of sites
	int<lower=0> nrep;									// Number of segments per site (always 5)
	int<lower=0> npred;									// Number of predictors
	
	matrix[nsite, nspec] X;								// Observed presence/absence matrix
	matrix[nsite, npred] preds;							// Predictor matrix 
}

parameters {
	matrix[nsite, nspec] Z;								// Latent occupancy matrix (what we are really interested in)
	
	real<lower=0> psi_mean;								// Hyperparam. of community-level occupancy
	real<lower=0> theta_mean;							// Hyperparam. of community-level detection
	vector[npred] mu_alpha;								// Hyperparam. of community-level habitat (alpha) and sampling (beta) covariates
	
	real<lower=0> tau1;
	real<lower=0> tau2;
	vector<lower=0>[npred] tau_alpha;
	real rho;
	
	vector[nspec] u;									// Occupancy covariate for each species
	vector[nspec] v;									// Detection covariate for each species
	matrix[npred, nspec] alpha;							// Habitat covariate(s) for each species
	
	matrix[nspec, nsite] psi;
}

transformed parameters {
	real a;
	real b;
	real var_v;
	real sigma1;
	real sigma2;
	vector[nspec] mu_v;
	
	a = log(psi_mean) - log(1 - psi_mean);
	b = log(theta_mean) - log(1 - theta_mean);
	var_v = tau2 / (1.0 - rho^2);
	sigma1 = 1 / sqrt(tau1);
	sigma2 = 1 / sqrt(tau2);
	mu_v = b + (rho * sigma2 / sigma1) * (u - a);
	
	// GENERATE psi HERE
}

model {
	// Likelihood
	for (i in 1:nsite) {
		Z[, i] ~ bernoulli(psi[, i])
	}
	
	
	// Priors
	psi_mean ~ uniform(0.001, 0.99);
	theta_mean ~ uniform(0.001, 0.99);
	tau1 ~ gamma(10, 1);
	tau2 ~ gamma(10, 1);
	tau_alpha ~ gamma(10, 1);
	rho ~ uniform(-0.99, 0.99);
	u ~ dnorm(a, tau1);
	v ~ dnorm(mu_v, var_v);
	alpha ~ dnorm(mu_alpha, tau_alpha); // have to change to a multivariate normal!
}

generated quantities {
	real occ_sp[nspec];								// Mean occupancy of species across all sites
}