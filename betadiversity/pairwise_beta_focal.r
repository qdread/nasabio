# Function that takes community matrix and calculates pairwise distances between the focal plot (row 1) and all other rows, and returns the mean of each.

pairwise_beta_focal <- function(m, td=TRUE, pd=TRUE, fd=TRUE, abundance=TRUE, pddist = NULL, fddist = NULL) {
	
	# Initialize data structures
	beta_td_pairwise_pa <- beta_td_shannon <- beta_td_pairwise <- rep(NA, nrow(m) - 1)
	beta_pd_pairwise_pa <- beta_pd_nt_pa <- beta_pd_pairwise <- beta_pd_nt <- rep(NA, nrow(m) - 1)
	beta_fd_pairwise_pa <- beta_fd_nt_pa <- beta_fd_pairwise <- beta_fd_nt <- rep(NA, nrow(m) - 1)
	
	# Calculate pairwise distances between focal community and other communities for each indicated type of diversity.
	for (i in 2:nrow(m)) {
		mat_i <- m[c(1,i), ]
		if (td) {
			beta_td_pairwise_pa[i] <- vegdist(x = mat_i, binary = TRUE, method = 'jaccard')
		}
		if (td & abundance) {
			beta_td_shannon[i] <- d(abundances = mat_i, lev = 'beta', wts = FALSE, q = 1)
			beta_td_pairwise[i] <- vegdist(x = mat_i, binary = FALSE, method = 'jaccard')
		}
		if (pd) {  
			beta_pd_pairwise_pa[i] <- comdist(comm = mat_i, dis = pddist, abundance.weighted = FALSE)
			beta_pd_nt_pa[i] <- comdistnt(comm = mat_i, dis = pddist, abundance.weighted = FALSE)
		}
		if (pd & abundance) {
			beta_pd_pairwise[i] <- comdist(comm = mat_i, dis = pddist, abundance.weighted = TRUE)
			beta_pd_nt[i] <- comdistnt(comm = mat_i, dis = pddist, abundance.weighted = TRUE)
		}
		if (fd) {
			beta_fd_pairwise_pa[i] <- comdist(comm = mat_i, dis = fddist, abundance.weighted = FALSE)
			beta_fd_nt_pa[i] <- comdistnt(comm = mat_i, dis = fddist, abundance.weighted = FALSE)
		}
		if (fd & abundance) {
			beta_fd_pairwise[i] <- comdist(comm = mat_i, dis = fddist, abundance.weighted = TRUE)
			beta_fd_nt[i] <- comdistnt(comm = mat_i, dis = fddist, abundance.weighted = TRUE)
		}
	}
	
	# Calculate means of each of the indicated diversity types.
	res <- c()
	
	if (td) {
		res <- c(res, beta_td_pairwise_pa = mean(beta_td_pairwise_pa, na.rm = TRUE))
	}
	if (td & abundance) {
		res <- c(res, beta_td_shannon = mean(beta_td_shannon, na.rm = TRUE),
					  beta_td_pairwise_pa = mean(beta_td_pairwise_pa, na.rm = TRUE))
	}
	if (pd) {
		res <- c(res, beta_pd_pairwise_pa = mean(beta_pd_pairwise_pa, na.rm = TRUE),
					  beta_pd_nt_pa = mean(beta_pd_nt_pa, na.rm = TRUE))
	}
	if (pd & abundance) {
		res <- c(res, beta_pd_pairwise = mean(beta_pd_pairwise, na.rm = TRUE),
					  beta_pd_nt = mean(beta_pd_nt, na.rm = TRUE))
	}
	if (fd) {
		res <- c(res, beta_fd_pairwise_pa = mean(beta_fd_pairwise_pa, na.rm = TRUE),
					  beta_fd_nt_pa = mean(beta_fd_nt_pa, na.rm = TRUE))
	}
	if (fd & abundance) {
		res <- c(res, beta_fd_pairwise = mean(beta_fd_pairwise, na.rm = TRUE),
					  beta_fd_nt = mean(beta_fd_nt, na.rm = TRUE))
	}
	return(res)
}