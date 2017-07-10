// Model predicting trait values from phylogenetic, environmental (and possibly spatial) information
// To be used for imputation
// First version will have spatial coordinates as fixed effects, not random

data {
	int m; 			// Number of traits
	int n;			// Number of species
	int N;			// Number of individuals across all species
	int p;			// Number of predictor variables + 1 (includes intercept)
	
	vector[N] Y;		// Trait vector, ordered by 1...m traits within each individual within each species
	matrix[m*N, p] X;	// Predictor matrix
	matrix[m*N, m*n] Z;	// Design matrix to map individuals to the proper species
	matrix[n, n] R;		// Phylogenetic variance-covariance matrix

}

transformed data {
	matrix[m*N, m*N] ImN;	// identity matrix dimension m*N
	matrix[m, m] Im;		// identity matrix dimension m
	vector[m*n] zero_mn;	// vector of zeroes m*n long
	vector[m*N] zero_mN;	// vector of zeroes m*N long
	
	zero_mn = rep_vector(0.0, m*n);
	zero_mN = rep_vector(0.0, m*N);
	Im = diag_matrix(rep_vector(1.0, m));
	ImN = diag_matrix(rep_vector(1.0, m * N));
}

parameters {
	vector[p] beta;				// Vector of coefficients on environmental predictors (fixed effects)
	vector[m*n] alpha;			// Vector of random effects
	
	matrix[m, m] Sigma;		// Trait variance-covariance matrix, to be estimated.
	real sigma;				// Hyperparameter on error term -- should this actually be a vector? Right now all individuals will have the same error.
}

transformed parameters {
	
	matrix[m*N, m*N] epsilon;	// Error term. Each individual gets its own error for now. Should have a diagonal of all sigma^2.
	matrix[m*n, m*n] lambda;
	
	epsilon = ImN * sigma^2;
	
	// Manually calculate Kronecker product of the trait covariance and the phylogenetic covariance matrices.
	for (i in 1:m)
		for (j in 1:m)
			for (k in 1:n)
				for (l in 1:n)
					lambda[n*(i-1)+k, n*(j-1)+l] = Sigma[i, j] * R[k, l];
}

model {
	// Likelihood
	Y ~ multi_normal(X * beta + Z * alpha, epsilon); 
	// Priors
	alpha ~ multi_normal(zero_mn, lambda);	// lambda is Kronecker product of trait and phylo covariance matrices (see below)
	beta ~ normal(0, 0.0001);
	Sigma ~ inv_wishart(m + 1, 0.01 * Im);
	sigma ~ uniform(0, 100);
}
