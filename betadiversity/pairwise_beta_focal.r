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
# Modification 19 May: standardize so that pd and fd are done with the same two distance-based metrics: pairwise and nearest-neighbor.

diversity_3ways <- function(m, flavor = 'alpha', dotd=TRUE, dopd=TRUE, dofd=TRUE, abundance=TRUE, pddist = NULL, fddist = NULL, nnull = 99, phylo_spp = NULL, func_problem_spp = NULL, combine = TRUE) {
	
  require(vegan)
  require(picante)
  
	# Get rid of species that aren't in phylogenetic and functional diversity.
	if (!is.null(phylo_spp)) m <- m[, dimnames(m)[[2]] %in% phylo_spp, drop = FALSE]
  if (!is.null(func_problem_spp)) mfunc <- m[, !dimnames(m)[[2]] %in% func_problem_spp, drop = FALSE] else mfunc <- m
	
	# Alpha diversity is done for each row separately.
	# Gamma diversity is for the sum of all the species down each column.
	if (flavor == 'gamma') {
		m <- t(apply(m, 2, sum))
		mfunc <- t(apply(mfunc, 2, sum))
	}
	
	# Declare variables to hold the data
	richness <- shannon <- evenness <- rep(NA, nrow(m))
	MPD_z <- MNTD_z <- rep(NA, nrow(m))
	MPD_pa_z <- MNTD_pa_z <- rep(NA, nrow(m))
	MPDfunc_z <- MNTDfunc_z <- rep(NA, nrow(m))
	MPDfunc_pa_z <- MNTDfunc_pa_z <- rep(NA, nrow(m))
	
	# Calculate the different flavors of diversity.
	if (dotd) {
		richness <- apply(m > 0, 1, sum)
	}
	if (dotd & abundance) {
		shannon <- diversity(m, index = 'shannon')
		evenness <- shannon/log(richness)
	}
	if (dopd) {
		MPD_pa <- ses.mpd(m, pddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
		MNTD_pa <- ses.mntd(m, pddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
		
		if (flavor == 'alpha') {
		  MPD_pa_z <- MPD_pa$mpd.obs.z
		  MNTD_pa_z <- MNTD_pa$mntd.obs.z
		}
		if (flavor == 'gamma') {
		  MPD_pa_z <- (MPD_pa$mpd.obs[1] - mean(MPD_pa$mpd.rand.mean, na.rm=T))/sd(MPD_pa$mpd.rand.mean, na.rm=T)
		  MNTD_pa_z <- (MNTD_pa$mntd.obs[1] - mean(MNTD_pa$mntd.rand.mean, na.rm=T))/sd(MNTD_pa$mntd.rand.mean, na.rm=T)
		}

	}
	if (dopd & abundance) {
		MPD <- ses.mpd(m, pddist, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = nnull)
		MNTD <- ses.mntd(m, pddist, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = nnull)
		if (flavor == 'alpha') {
		  MPD_z <- MPD$mpd.obs.z
		  MNTD_z <- MNTD$mntd.obs.z
		}
		if (flavor == 'gamma') {
		  MPD_z <- (MPD$mpd.obs[1] - mean(MPD$mpd.rand.mean, na.rm=T))/sd(MPD$mpd.rand.mean, na.rm=T)
		  MNTD_z <- (MNTD$mntd.obs[1] - mean(MNTD$mntd.rand.mean, na.rm=T))/sd(MNTD$mntd.rand.mean, na.rm=T)
		}
	}
	if (dofd) {
	  MPDfunc_pa <- ses.mpd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
	  MNTDfunc_pa <- ses.mntd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
	  
	  if (flavor == 'alpha') {
	    MPDfunc_pa_z <- MPDfunc_pa$mpd.obs.z
	    MNTDfunc_pa_z <- MNTDfunc_pa$mntd.obs.z
	  }
	  if (flavor == 'gamma') {
	    MPDfunc_pa_z <- (MPDfunc_pa$mpd.obs[1] - mean(MPDfunc_pa$mpd.rand.mean, na.rm=T))/sd(MPDfunc_pa$mpd.rand.mean, na.rm=T)
	    mntdrand <- MNTDfunc_pa$mntd.rand.mean[is.finite(MNTDfunc_pa$mntd.rand.mean)]
	    MNTDfunc_pa_z <- (MNTDfunc_pa$mntd.obs[1] - mean(mntdrand))/sd(mntdrand)
	  }
	  
	}
	if (dofd & abundance) {
	  MPDfunc <- ses.mpd(mfunc, pddist, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = nnull)
	  MNTDfunc <- ses.mntd(mfunc, pddist, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = nnull)
	  if (flavor == 'alpha') {
	    MPDfunc_z <- MPDfunc$mpd.obs.z
	    MNTDfunc_z <- MNTDfunc$mntd.obs.z
	  }
	  if (flavor == 'gamma') {
	    MPDfunc_z <- (MPDfunc$mpd.obs[1] - mean(MPDfunc$mpd.rand.mean, na.rm=T))/sd(MPDfunc$mpd.rand.mean, na.rm=T)
	    mntdrand <- MNTDfunc$mntd.rand.mean[is.finite(MNTDfunc$mntd.rand.mean)]
	    MNTDfunc_z <- (MNTDfunc$mntd.obs[1] - mean(mntdrand))/sd(mntdrand)
	  }
	}
	
	
	# Concatenate into a vector and return.
	if(combine) return(c(richness = median(richness, na.rm=T), shannon = median(shannon, na.rm=T), evenness = median(evenness, na.rm=T),
	  MPD_pa_z = median(MPD_pa_z, na.rm=T), MNTD_pa_z = median(MNTD_pa_z, na.rm=T),
	  MPD_z = median(MPD_z, na.rm=T), MNTD_z = median(MNTD_z, na.rm=T),
	  MPDfunc_pa_z = median(MPDfunc_pa_z, na.rm=T), MNTDfunc_pa_z = median(MNTDfunc_pa_z, na.rm=T),
	  MPDfunc_z = median(MPDfunc_z, na.rm=T), MNTDfunc_z = median(MNTDfunc_z, na.rm=T)))
	if(!combine) return(data.frame(richness, shannon, evenness, MPD_pa_z, MNTD_pa_z, MPD_z, MNTD_z, MPDfunc_pa_z, MNTDfunc_pa_z, MPDfunc_z, MNTDfunc_z))
	

}