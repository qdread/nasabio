// Jarzyna's occupancy model translated from JAGS to STAN
// QDR 17 July 2017
// Goal is to speed up the fitting process (vectorized instead of with loops)
// Also generalize so that we can include multiple predictors.

data {
	int<lower=0> nspec;									// Number of species
	int<lower=0> nsite;									// Number of sites
	int<lower=0> nrep;									// Number of segments per site (always 5)
	int<lower=0> npred;									// Number of predictors
	
	int<lower=0, upper=1> X[nsite, nrep, nspec];		// Observed presence/absence array
	matrix[nsite, npred] preds;							// Predictor matrix 
}

parameters {
	
	real<lower=0> psi_mean;								// Hyperparam. of community-level occupancy
	real<lower=0> theta_mean;							// Hyperparam. of community-level detection
	vector[npred] mu_alpha;								// Hyperparam. of community-level habitat (alpha) and sampling (beta) covariates
	
	real<lower=0> tau1;
	real<lower=0> tau2;
	cov_matrix[npred] sigma_alpha;						// Changed from tau
	real rho;
	
	vector[nspec] u;									// Occupancy covariate for each species
	vector[nspec] v;									// Detection covariate for each species
	matrix[npred, nspec] alpha;							// Habitat covariate(s) for each species
	
}

transformed parameters {
	real a;
	real b;
	real var_v;
	real<lower=0> sigma1;
	real<lower=0> sigma2;
	vector[nspec] mu_v;
	matrix<lower=0, upper=1>[nsite, nspec] psi;			// Occurrence probability matrix
	real<lower=0, upper=1> theta[nsite, nrep, nspec];	// Detection probability array

	a = log(psi_mean) - log(1 - psi_mean);
	b = log(theta_mean) - log(1 - theta_mean);
	var_v = tau2 / (1.0 - rho^2);
	sigma1 = 1 / sqrt(tau1);
	sigma2 = 1 / sqrt(tau2);
	mu_v = b + (rho * sigma2 / sigma1) * (u - a);
	
	for (i in 1:nspec) {
		for (j in 1:nsite) {
			
			psi[j, i] = inv_logit(u[i] + preds[i, ] * alpha[, i]);
			
			for (k in 1:nrep) {
				theta[j, k, i] = inv_logit(v[i]);
			}
		}
	}
}

model {
	// Likelihood
	for (i in 1:nspec) {
		for (j in 1:nsite) {
			for (k in 1:nrep) {
				X[j, k, i] ~ bernoulli(theta[j, k, i] * psi[j, i]); // last term is Z[j,i] from old model
			}
		}
	}

	// Priors
	psi_mean ~ uniform(0.001, 0.99);
	theta_mean ~ uniform(0.001, 0.99);
	tau1 ~ gamma(10, 1);
	tau2 ~ gamma(10, 1);
	mu_alpha ~ normal(0, 100);
	rho ~ uniform(-0.99, 0.99);
	u ~ normal(a, sigma1);
	v ~ normal(mu_v, var_v^-2);
	for (i in 1:nspec) {
		alpha[, i] ~ multi_normal(mu_alpha, sigma_alpha);
	}
}
