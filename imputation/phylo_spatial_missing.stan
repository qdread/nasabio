// Model predicting trait values from phylogenetic, environmental (and possibly spatial) information
// This version includes the possibility of missing values in Y.
// See Stan reference section 11.3: Sliced Missing Data
// First version will have spatial coordinates as fixed effects, not random

data {
	int m; 			// Number of traits
	int n;			// Number of species
	int N;			// Number of individuals across all species
	int p;			// Number of predictor variables + 1 (includes intercept)
	int N_obs;		// Number of values in Y that are observed (not missing)
	int N_mis;		// Number of values in Y that are missing
					// ***NOTE: n_obs + n_mis must be m*N !!!
	
	vector[N_obs] Y_obs;	// Trait vector, ordered by 1...m traits within each individual within each species. Shortened by the number of missing values.
	matrix[m*N, m*p] X;		// Predictor matrix
	matrix[m*N, m*n] Z;		// Design matrix to map individuals to the proper species
	cov_matrix[n] R;		// Phylogenetic variance-covariance matrix
	
	int<lower = 1, upper = N_obs + N_mis> index_obs[N_obs];		// Indexes of observed values
	int<lower = 1, upper = N_obs + N_mis> index_mis[N_mis];		// Indexes of missing values
}

transformed data {
	matrix[m*N, m*N] ImN;	// identity matrix dimension m*N
	matrix[m, m] Im;		// identity matrix dimension m
	vector[m*n] zero_mn;	// vector of zeroes m*n long
	vector[m*N] zero_mN;	// vector of zeroes m*N long
	
	zero_mn = rep_vector(0.0, m*n);
	Im = diag_matrix(rep_vector(1.0, m));
}

parameters {
	vector[m*p] beta;			// Vector of coefficients on environmental predictors (fixed effects)
	vector[m*n] alpha;			// Vector of random effects
	
	cov_matrix[m] Sigma;		// Trait variance-covariance matrix, to be estimated.
	vector[m] sigma;			// Hyperparameter on error term--different error for each trait.
	
	vector[N_mis] Y_mis;		// Missing values as parameters to be estimated.
}

transformed parameters {
	
	vector[m*N] Y;				// Y will be stitched back together from Y_obs and Y_mis.
	
	cov_matrix[m*N] epsilon;	// Error term. Each individual gets its own error for now. Should have a diagonal of all sigma^2s.
	vector[m*N] epsilon_diag;	// Diagonal of epsilon matrix.	
	vector[m] sigma2;			

	cov_matrix[m*n] lambda;		// Kronecker product of Sigma (trait vcov matrix) and R (phylogenetic vcov matrix)
	
	for (i in 1:m)
		sigma2[i] = sigma[i]^2;
	
	for (i in 1:m)
		epsilon_diag[(1 + (N*(i-1))):(N + (N*(i-1)))] = rep_vector(sigma2[i], N);

	epsilon = diag_matrix(epsilon_diag);
	
	// Manually calculate Kronecker product of the trait covariance and the phylogenetic covariance matrices.
	for (i in 1:m)
		for (j in 1:m)
			for (k in 1:n)
				for (l in 1:n)
					lambda[n*(i-1)+k, n*(j-1)+l] = Sigma[i, j] * R[k, l];
				
	// Re-assemble Y including both missing and observed values.
	Y[index_obs] = Y_obs;
	Y[index_mis] = Y_mis;
}

model {
	// Likelihood
	Y ~ multi_normal(X * beta + Z * alpha, epsilon);
	// Priors
	alpha ~ multi_normal(zero_mn, lambda);	
	beta ~ normal(0, 10);
	Sigma ~ inv_wishart(m + 1, 0.01 * Im);
	sigma ~ uniform(0, 100);
}
