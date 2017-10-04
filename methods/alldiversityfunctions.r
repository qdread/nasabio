# Functions for calculating diversity indices
# alpha, beta, and gamma x taxonomic, functional, phylogenetic (9 combinations)
# Function diversity_3ways() takes a set of communities as input. It can calculate td, pd, and fd for each of the communities separately if "alpha" is chose
# If "gamma" is chosen, it aggregates the communities and calculates td, pd, and fd
# Later, the individual alpha values can be aggregated by radius. Gamma is already aggregated.

# Function singlepair_beta() takes two communities and calculates pairwise beta diversity between them, for td, pd, and fd
# For both pd and fd this is distance-based.
# This is used to generate a lookup table of pairwise measurements that can be averaged at a desired radius
# Function beta_part() takes a set of communities and calculates the multisite Sorensen based beta diversity for all of them. Unfortunately this is only possible for td and pd, not fd.
# As an added bonus, it also calculates the Baselga partition of each type of beta-diversity into nestedness and turnover components.

# Last updated by QDR on 04 Oct 2017
# NASABIOXGEO Project

#####################################################################
# diversity_3ways(): function to do all 3 flavors of diversity for alpha or gamma given an input community matrix.
# It will return many indices, one for each community if flavor='alpha', one for the entire set if flavor='gamma'

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


########################################################################
# singlepair_beta(): pairwise beta for a single pair of sites.
# Also does the null model.

singlepair_beta <- function(p1, p2, td=TRUE, pd=TRUE, fd=TRUE, abundance=TRUE, pddist = NULL, fddist = NULL, nnull = 99, phylo_spp = NULL, func_problem_spp = NULL) {
	
	m <- rbind(p1, p2)
	
	# Get rid of species that aren't in phylogenetic and functional diversity.
	if (!is.null(phylo_spp)) m <- m[, dimnames(m)[[2]] %in% phylo_spp, drop = FALSE]
    if (!is.null(func_problem_spp)) mfunc <- m[, !dimnames(m)[[2]] %in% func_problem_spp, drop = FALSE] else mfunc <- m
	
	# Declare variables to hold the data
	beta_td_pairwise_pa <- beta_td_sorensen_pa <- beta_td_shannon <- beta_td_pairwise <- beta_td_sorensen <- NA
	beta_pd_pairwise_pa <- beta_pd_nt_pa <- beta_pd_pairwise <- beta_pd_nt <- NA
	beta_fd_pairwise_pa <- beta_fd_nt_pa <- beta_fd_pairwise <- beta_fd_nt <- NA
	
	# Calculate pairwise distances between focal community and other communities for each indicated type of diversity.

	if (td) {
		beta_td_pairwise_pa <- vegdist(x = m, binary = TRUE, method = 'jaccard')[1]
		beta_td_sorensen_pa <- vegdist(x = m, binary = TRUE, method = 'bray')[1]
	}
	if (td & abundance) {
		beta_td_shannon <- d(abundances = m, lev = 'beta', wts = FALSE, q = 1)
		beta_td_pairwise <- vegdist(x = m, binary = FALSE, method = 'jaccard')[1]				
		beta_td_sorensen <- vegdist(x = m, binary = FALSE, method = 'bray')[1]
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
	c(beta_td_pairwise_pa=beta_td_pairwise_pa, beta_td_sorensen_pa=beta_td_sorensen_pa,
	  beta_td_pairwise=beta_td_pairwise, beta_td_sorensen=beta_td_sorensen, beta_td_shannon=beta_td_shannon,
	  beta_pd_pairwise_pa=beta_pd_pairwise_pa, beta_pd_pairwise_pa_z=beta_pd_pairwise_pa_z,
	  beta_pd_nt_pa=beta_pd_nt_pa, beta_pd_nt_pa_z=beta_pd_nt_pa_z,
	  beta_pd_pairwise=beta_pd_pairwise, beta_pd_pairwise_z=beta_pd_pairwise_z,
	  beta_pd_nt=beta_pd_nt, beta_pd_nt_pa_z=beta_pd_nt_z,
	  beta_fd_pairwise_pa=beta_fd_pairwise_pa, beta_fd_pairwise_pa_z=beta_fd_pairwise_pa_z,
	  beta_fd_nt_pa=beta_fd_nt_pa, beta_fd_nt_pa_z=beta_fd_nt_pa_z,
	  beta_fd_pairwise=beta_fd_pairwise, beta_fd_pairwise_z=beta_fd_pairwise_z,
	  beta_fd_nt=beta_fd_nt, beta_fd_nt_pa_z=beta_fd_nt_z)

}

######################################################################
# beta_part(): Takes a set of communities and does multi-site beta-diversity
# It also does the Baselga partition for both.

beta_part <- function(m, abundance=TRUE, pairwise=FALSE, index_family='sorensen', TD=TRUE, FD=FALSE, PD=FALSE, trait_mat = NULL, phylo_tree = NULL) {
  
  require(betapart)
  m_bin <- m
  m_bin[m_bin > 0] <- 1
  
  index_family_abund <- ifelse(index_family == 'sorensen', 'bray', 'ruzicka')
  
  # Do precalculation.
  core_presence  <- betapart.core(m_bin)
  
  # Calculate metrics.
  beta_presence  <- unlist(beta.multi(core_presence, index.family = index_family))
  
  if (abundance) {
    core_abundance <- betapart.core.abund(m)
    beta_abundance <- unlist(beta.multi.abund(core_abundance, index.family = index_family_abund))
  }
  
  if (FD) {
    trait_mat <- trait_mat[dimnames(trait_mat)[[1]] %in% dimnames(m_bin)[[2]], ]
    core_func <- functional.betapart.core(m_bin, traits = as.matrix(trait_mat))
    beta_func <- unlist(functional.beta.multi(core_func, index.family = index_family))
  }
  
  if (PD) {
    core_phy <- phylo.betapart.core(m_bin, tree = phylo_tree)
    beta_phy <- unlist(phylo.beta.multi(core_phy, index.family = index_family))
  }
  
  # Calculate pairwise metrics if needed.
  if (pairwise) {
    beta_presence_pair  <- beta.pair(core_presence, index.family = index_family)
    if (abundance) {
      beta_abundance_pair <- beta.pair.abund(core_abundance, index.family = index_family_abund)
    }
    if (FD) {
      beta_func_pair <- functional.beta.pair(core_func, traits = trait_mat, index.family = index_family)
    }
    if (PD) {
      beta_phylo_pair <- phylo.beta.pair(core_phy, tree = phylo_tree, index.family = index_family)
    }
  }

  # Combine and return results.
  res <- data.frame(index = rep(index_family, 3),
                    diversity = 'taxonomic',
                    partition = c('replacement','nestedness','total'),
                    abundance = FALSE,
                    beta = beta_presence)
  if (abundance) {
    res_abund <- data.frame(index = rep(index_family_abund, 3),
                            diversity = 'taxonomic',
                            partition = c('replacement','nestedness','total'),
                            abundance = TRUE,
                            beta = beta_abundance)
    res <- rbind(res, res_abund)
  }
  if (FD) {
    res_func <- data.frame(index = rep(index_family, 3),
                           diversity = 'functional',
                           partition = c('replacement','nestedness','total'),
                           abundance = FALSE,
                           beta = beta_func)
    res <- rbind(res, res_func)
  }
  if (PD) {
    res_phy <- data.frame(index = rep(index_family, 3),
                          diversity = 'phylogenetic',
                          partition = c('replacement','nestedness','total'),
                          abundance = FALSE,
                          beta = beta_phy)
    res <- rbind(res, res_phy)
  }
  
  if (!pairwise) return(res)
  res <- list(res, beta_presence_pair)
  if (abundance) res[[length(res) + 1]] <- beta_abundance_pair
  if (FD) res[[length(res) + 1]] <- beta_func_pair
  if (PD) res[[length(res) + 1]] <- beta_phy_pair
  return(res)
}