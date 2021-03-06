// Jarzyna's occupancy model translated from JAGS to STAN
// QDR 20 July 2017
// Goal is to speed up the fitting process (vectorized instead of with loops)
// This one only includes random effects and has no fixed predictors.
// This version has less parameters to be monitored, to save space.

data {
	int<lower=0> nspec;									// Number of species
	int<lower=0> nsite;									// Number of sites
	int<lower=0> nrep;									// Number of segments per site (always 5)
	
	int<lower=0, upper=1> X[nsite, nrep, nspec];		// Observed presence/absence array
}

parameters {
	
	real<lower=0> psi_mean;								// Hyperparam. of community-level occupancy
	real<lower=0> theta_mean;							// Hyperparam. of community-level detection
	
	real<lower=0> tau1;
	real<lower=0> tau2;
	real<lower=-1, upper=1> rho;
	
	vector[nspec] u;									// Occupancy covariate for each species
	vector[nspec] v;									// Detection covariate for each species
	
}

transformed parameters {
	real a;
	real b;
	real<lower=0> sigma_v;
	real<lower=0> sigma1;
	real<lower=0> sigma2;
	vector[nspec] mu_v;

	a = log(psi_mean) - log(1 - psi_mean);
	b = log(theta_mean) - log(1 - theta_mean);
	sigma_v = (1.0 - rho^2) / tau2;
	sigma1 = 1 / sqrt(tau1);
	sigma2 = 1 / sqrt(tau2);
	mu_v = b + (rho * sigma2 / sigma1) * (u - a);
}

model {
	// Likelihood
	for (i in 1:nspec) {
		for (j in 1:nsite) {
			for (k in 1:nrep) {
				// P(bird is seen) = P(bird is there) * P(bird is seen | bird is there)
				X[j, k, i] ~ bernoulli(inv_logit(v[i]) * inv_logit(u[i])); 
			}
		}
	}

	// Priors
	psi_mean ~ uniform(0.001, 0.99);
	theta_mean ~ uniform(0.001, 0.99);
	tau1 ~ gamma(10, 1);
	tau2 ~ gamma(10, 1);
	rho ~ uniform(-0.99, 0.99);
	u ~ normal(a, sigma1);
	v ~ normal(mu_v, sigma_v);
}
