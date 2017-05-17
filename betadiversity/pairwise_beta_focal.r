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


#################################################################
# Added 17 May: pairwise beta for a single pair.
# Also does the null model.

singlepair_beta <- function(p1, p2, td=TRUE, pd=TRUE, fd=TRUE, abundance=TRUE, pddist = NULL, fddist = NULL, nnull = 99, phylo_spp = NULL, func_problem_spp = NULL) {
	
	m <- rbind(p1, p2)
	
	# Get rid of species that aren't in phylogenetic and functional diversity.
	if (!is.null(phylo_spp)) m <- m[, dimnames(m)[[2]] %in% phylo_spp, drop = FALSE]
    if (!is.null(func_problem_spp)) mfunc <- m[, !dimnames(m)[[2]] %in% func_problem_spp, drop = FALSE] else mfunc <- m
	
	# Declare variables to hold the data
	beta_td_pairwise_pa <- beta_td_shannon <- beta_td_pairwise <- NA
	beta_pd_pairwise_pa <- beta_pd_nt_pa <- beta_pd_pairwise <- beta_pd_nt <- NA
	beta_fd_pairwise_pa <- beta_fd_nt_pa <- beta_fd_pairwise <- beta_fd_nt <- NA
	
	# Calculate pairwise distances between focal community and other communities for each indicated type of diversity.

	if (td) {
		beta_td_pairwise_pa <- vegdist(x = m, binary = TRUE, method = 'jaccard')[1]
	}
	if (td & abundance) {
		beta_td_shannon <- d(abundances = m, lev = 'beta', wts = FALSE, q = 1)
		beta_td_pairwise <- vegdist(x = m, binary = FALSE, method = 'jaccard')[1]
	}
	if (pd) {  
		beta_pd_pairwise_pa <- comdist(comm = m, dis = pddist, abundance.weighted = FALSE)[1]
		beta_pd_nt_pa <- comdistnt(comm = m, dis = pddist, abundance.weighted = FALSE)[1]
	}
	if (pd & abundance) {
		beta_pd_pairwise <- comdist(comm = m, dis = pddist, abundance.weighted = TRUE)[1]
		beta_pd_nt <- comdistnt(comm = m, dis = pddist, abundance.weighted = TRUE)[1]
	}
	if (fd & ncol(mfunc) > 1) {
		beta_fd_pairwise_pa <- comdist(comm = mfunc, dis = fddist, abundance.weighted = FALSE)[1]
		beta_fd_nt_pa <- comdistnt(comm = mfunc, dis = fddist, abundance.weighted = FALSE)[1]
	}
	if (fd & abundance & ncol(mfunc) > 1) {
		beta_fd_pairwise <- comdist(comm = mfunc, dis = fddist, abundance.weighted = TRUE)[1]
		beta_fd_nt <- comdistnt(comm = mfunc, dis = fddist, abundance.weighted = TRUE)[1]
	}
		
	# Get null distribution of pd and fd beta-diversity metrics and calculate the z-scores for each.
	phypairwise_null <- phypairwise_pa_null <- phynt_null <- phynt_pa_null <- rep(NA, nnull)
	funcpairwise_null <- funcpairwise_pa_null <- funcnt_null <- funcnt_pa_null <- rep(NA, nnull)
		for (sim in 1:nnull) {
			if (pd) {
						
				pdnullidx <- sample(1:nrow(pddist))
				pddistnull <- pddist
				dimnames(pddistnull)[[1]] <- dimnames(pddistnull)[[1]][pdnullidx]
				dimnames(pddistnull)[[2]] <- dimnames(pddistnull)[[2]][pdnullidx]
							
				phypairwise_pa_null[sim] <- mean(comdist(comm = m, dis = pddistnull, abundance.weighted = FALSE))
				phynt_pa_null[sim] <- mean(comdistnt(comm = m, dis = pddistnull, abundance.weighted = FALSE))
				
				if (abundance) {
					phypairwise_null[sim] <- mean(comdist(comm = m, dis = pddistnull, abundance.weighted = TRUE))
					phynt_null[sim] <- mean(comdistnt(comm = m, dis = pddistnull, abundance.weighted = TRUE))
				}
						
            }
			
			if (fd & ncol(mfunc) > 1) {
			 
				fdnullidx <- sample(1:nrow(fddist))
				fddistnull <- fddist
				dimnames(fddistnull)[[1]] <- dimnames(fddistnull)[[1]][fdnullidx]
				dimnames(fddistnull)[[2]] <- dimnames(fddistnull)[[2]][fdnullidx]
				
				funcpairwise_pa_null[sim] <- mean(comdist(comm = mfunc, dis = fddistnull, abundance.weighted = FALSE))
				funcnt_pa_null[sim] <- mean(comdistnt(comm = mfunc, dis = fddistnull, abundance.weighted = FALSE))
				
				if (abundance) {
					funcpairwise_null[sim] <- mean(comdist(comm = mfunc, dis = fddistnull, abundance.weighted = TRUE))
					funcnt_null[sim] <- mean(comdistnt(comm = mfunc, dis = fddistnull, abundance.weighted = TRUE))
				}
						
			}
		}
	
	# Calculation of z-scores.
	beta_pd_pairwise_z <- (beta_pd_pairwise - mean(phypairwise_null, na.rm=T))/sd(phypairwise_null, na.rm=T)
    beta_pd_pairwise_pa_z <- (beta_pd_pairwise_pa - mean(phypairwise_pa_null, na.rm=T))/sd(phypairwise_pa_null, na.rm=T)
    beta_pd_nt_z <- (beta_pd_nt - mean(phynt_null, na.rm=T))/sd(phynt_null, na.rm=T)
    beta_pd_nt_pa_z <- (beta_pd_nt_pa - mean(phynt_pa_null, na.rm=T))/sd(phynt_pa_null, na.rm=T)
	beta_fd_pairwise_z <- (beta_fd_pairwise - mean(funcpairwise_null, na.rm=T))/sd(funcpairwise_null, na.rm=T)
    beta_fd_pairwise_pa_z <- (beta_fd_pairwise_pa - mean(funcpairwise_pa_null, na.rm=T))/sd(funcpairwise_pa_null, na.rm=T)
    beta_fd_nt_z <- (beta_fd_nt - mean(funcnt_null, na.rm=T))/sd(funcnt_null, na.rm=T)
    beta_fd_nt_pa_z <- (beta_fd_nt_pa - mean(funcnt_pa_null, na.rm=T))/sd(funcnt_pa_null, na.rm=T)
	
	# Concatenate results into a vector and return them.
	c(beta_td_pairwise_pa=beta_td_pairwise_pa, beta_td_pairwise=beta_td_pairwise, beta_td_shannon=beta_td_shannon, 
	  beta_pd_pairwise_pa=beta_pd_pairwise_pa, beta_pd_pairwise_pa_z=beta_pd_pairwise_pa_z,
	  beta_pd_nt_pa=beta_pd_nt_pa, beta_pd_nt_pa_z=beta_pd_nt_pa_z,
	  beta_pd_pairwise=beta_pd_pairwise, beta_pd_pairwise_z=beta_pd_pairwise_z,
	  beta_pd_nt=beta_pd_nt, beta_pd_nt_pa_z=beta_pd_nt_z,
	  beta_fd_pairwise_pa=beta_fd_pairwise_pa, beta_fd_pairwise_pa_z=beta_fd_pairwise_pa_z,
	  beta_fd_nt_pa=beta_fd_nt_pa, beta_fd_nt_pa_z=beta_fd_nt_pa_z,
	  beta_fd_pairwise=beta_fd_pairwise, beta_fd_pairwise_z=beta_fd_pairwise_z,
	  beta_fd_nt=beta_fd_nt, beta_fd_nt_pa_z=beta_fd_nt_z)

}


#####################################################################
# Added 17 May: function to do all 3 flavors of diversity for alpha or gamma given an input matrix. The matrix is all plots within a given radius of the focal point.

diversity_3ways <- function(m, flavor = 'alpha', td=TRUE, pd=TRUE, fd=TRUE, abundance=TRUE, pddist = NULL, fddist = NULL, nnull = 99, phylo_spp = NULL, func_problem_spp = NULL, phy = NULL) {
	
	# Get rid of species that aren't in phylogenetic and functional diversity.
	if (!is.null(phylo_spp)) m <- m[, dimnames(m)[[2]] %in% phylo_spp, drop = FALSE]
    if (!is.null(func_problem_spp)) mfunc <- m[, !dimnames(m)[[2]] %in% func_problem_spp, drop = FALSE] else mfunc <- m
	
	# Alpha diversity is done for each row separately.
	# Gamma diversity is for the sum of all the species down each column.
	if (flavor == 'gamma') {
		m <- matrix(apply(m, 2, sum), nrow=1)
		mfunc <- matrix(apply(mfunc, 2, sum), nrow=1)
	}
	
	# Declare variables to hold the data
	richness <- shannon <- evenness <- rep(NA, nrow(m))
	PD <- MPD_z <- MNTD_z <- rep(NA, nrow(m))
	PD_pa <- MPD_pa_z <- MNTD_pa_z <- rep(NA, nrow(m))
	FRic <- FDis <- FDiv <- FEve <- rep(NA, nrow(m))
	FRic_pa <- FDis_pa <- FDiv_pa <- FEve_pa <- rep(NA, nrow(m))
		
	# Calculate the different flavors of diversity.
	if (td) {
		richness <- apply(m > 0, 1, sum)
	}
	if (td & abundance) {
		shannon <- diversity(m, index = 'shannon')
		evenness <- shannon/log(richness)
	}
	if (pd) {
		PD_pa <- pd(m > 0, phy, include.root = TRUE)
		MPD_pa_z <- ses.mpd(m, pddist, null.model = 'independentswap', abundance.weighted = FALSE, runs = nnull, iterations = 1000)
		MNTD_pa_z <- ses.mntd(m, pddist, null.model = 'independentswap', abundance.weighted = FALSE, runs = nnull, iterations = 1000)
	}
	if (pd & abundance) {
		PD <- pd(m, phy, include.root = TRUE)
		MPD_z <- ses.mpd(m, pddist, null.model = 'independentswap', abundance.weighted = TRUE, runs = nnull, iterations = 1000)
		MNTD_z <- ses.mntd(m, pddist, null.model = 'independentswap', abundance.weighted = TRUE, runs = nnull, iterations = 1000)
	}
	if (fd) {
		zerorows <- apply(mfunc, 1, sum) > 0
		fd_plot_pa <- dbFD(x = , a = mfunc[!zerorows, ], w.abun = FALSE, corr = 'cailliez')

		# Output of FD should have NAs for the zero-abundance communities
		FRic_pa[!zerorows] <- fd_plot_pa$FRic
		FEve_pa[!zerorows] <- fd_plot_pa$FEve
		FDiv_pa[!zerorows] <- fd_plot_pa$FDiv
		FDis_pa[!zerorows] <- fd_plot_pa$FDis
	}
	if (fd & abundance) {
		fd_plot <- dbFD(x = , a = mfunc[!zerorows, ], w.abun = TRUE, corr = 'cailliez')

		# Output of FD should have NAs for the zero-abundance communities
		FRic[!zerorows] <- fd_plot$FRic
		FEve[!zerorows] <- fd_plot$FEve
		FDiv[!zerorows] <- fd_plot$FDiv
		FDis[!zerorows] <- fd_plot$FDis
	}
	
	
	# Concatenate into a vector and return.
	c(richness = median(richness, na.rm=T), shannon = median(shannon, na.rm=T), evenness = median(evenness, na.rm=T),
	  PD_pa = median(PD_pa, na.rm=T), MPD_pa_z = median(MPD_pa_z, na.rm=T), MNTD_pa_z = median(MNTD_pa_z, na.rm=T),
	  PD = median(PD, na.rm=T), MPD_z = median(MPD_z, na.rm=T), MNTD_z = median(MNTD_z, na.rm=T),
	  FRic = median(FRic, na.rm=T), FEve = median(FEve, na.rm=T), FDiv = median(FDiv, na.rm=T), FDis = median(FDis, na.rm=T),
	  FRic_pa = median(FRic_pa, na.rm=T), FEve_pa = median(FEve_pa, na.rm=T), FDiv_pa = median(FDiv_pa, na.rm=T), FDis_pa = median(FDis_pa, na.rm=T))
	

}