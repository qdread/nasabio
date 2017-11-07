radii <- c(1000,5000,7500,10000,20000,50000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

r <- radii[task]


load(paste0('/mnt/research/nasabio/data/fia/mats/mat_', r, '.r'))
fia_list <- list()

load('/mnt/research/nasabio/data/fia/fiaworkspace.r') 

library(vegan)
library(vegetarian)
source('~/code/fia/fixpicante.r')

nnull <- 99 # Reduce to save time

trydist <- as.matrix(trydist)

pb <- txtProgressBar(0, length(all_mats), style=3)

for (p in 1:length(all_mats)) {
	setTxtProgressBar(pb, p)
	mat_p <- all_mats[[p]]
	
	if(inherits(mat_p, 'matrix')) {
		  if(nrow(mat_p) > 1 & ncol(mat_p) > 1) {
			
			# Fix the species names to match the phylogeny, and get rid of the unknown species.
        
        mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% pnwphylo$tip.label, drop = FALSE]
        mat_p_noproblem <- mat_p[, !dimnames(mat_p)[[2]] %in% problemspp, drop = FALSE]
        
        # Calculate beta-diversity for that matrix.
        
        fia_shannonbetadiv <- d(abundances = mat_p, lev = 'beta', wts = FALSE, q = 1)
        fia_meanpairwisedissim <- mean(vegdist(x = mat_p, binary = FALSE, method = 'jaccard'))
        fia_meanpairwisedissim_pa <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        
			fia_phypairwise <- mean(comdist(comm = mat_p, dis = fiadist, abundance.weighted = TRUE))
          fia_phypairwise_pa <- mean(comdist(comm = mat_p, dis = fiadist, abundance.weighted = FALSE))
          fia_phynt <- mean(comdistnt(comm = mat_p, dis = fiadist, abundance.weighted = TRUE))
          fia_phynt_pa <- mean(comdistnt(comm = mat_p, dis = fiadist, abundance.weighted = FALSE))
		  if (ncol(mat_p_noproblem) > 1) {
			  fia_funcpairwise <- mean(comdist(comm = mat_p_noproblem, dis = trydist, abundance.weighted = TRUE))
			  fia_funcpairwise_pa <- mean(comdist(comm = mat_p_noproblem, dis = trydist, abundance.weighted = FALSE))
			  fia_funcnt <- mean(comdistnt(comm = mat_p_noproblem, dis = trydist, abundance.weighted = TRUE))
			  fia_funcnt_pa <- mean(comdistnt(comm = mat_p_noproblem, dis = trydist, abundance.weighted = FALSE))
		  }
		  else {
			fia_funcpairwise <- NA
			fia_funcpairwise_pa <- NA
			fia_funcnt <- NA
			fia_funcnt_pa <- NA
		  }
          # Null models by scrambling distance matrix
          phypairwise_null <- phypairwise_pa_null <- phynt_null <- phynt_pa_null <- rep(NA, nnull)
		  funcpairwise_null <- funcpairwise_pa_null <- funcnt_null <- funcnt_pa_null <- rep(NA, nnull)
          
          for (sim in 1:nnull) {
            nullidx <- sample(1:nrow(fiadist))
            fiadistnull <- fiadist
            dimnames(fiadistnull)[[1]] <- dimnames(fiadistnull)[[1]][nullidx]
            dimnames(fiadistnull)[[2]] <- dimnames(fiadistnull)[[2]][nullidx]
            
			trynullidx <- sample(1:nrow(trydist))
            trydistnull <- trydist
            dimnames(trydistnull)[[1]] <- dimnames(trydistnull)[[1]][trynullidx]
            dimnames(trydistnull)[[2]] <- dimnames(trydistnull)[[2]][trynullidx]
			
            phypairwise_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
            phypairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
            phynt_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
            phynt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
			
			if (ncol(mat_p_noproblem) > 1) {
			funcpairwise_null[sim] <- mean(comdist(comm = mat_p_noproblem, dis = trydistnull, abundance.weighted = TRUE))
            funcpairwise_pa_null[sim] <- mean(comdist(comm = mat_p_noproblem, dis = trydistnull, abundance.weighted = FALSE))
            funcnt_null[sim] <- mean(comdistnt(comm = mat_p_noproblem, dis = trydistnull, abundance.weighted = TRUE))
            funcnt_pa_null[sim] <- mean(comdistnt(comm = mat_p_noproblem, dis = trydistnull, abundance.weighted = FALSE))
			}
			
          }
          
          fia_phypairwise_z <- (fia_phypairwise - mean(phypairwise_null, na.rm=T))/sd(phypairwise_null, na.rm=T)
          fia_phypairwise_pa_z <- (fia_phypairwise_pa - mean(phypairwise_pa_null, na.rm=T))/sd(phypairwise_pa_null, na.rm=T)
          fia_phynt_z <- (fia_phynt - mean(phynt_null, na.rm=T))/sd(phynt_null, na.rm=T)
          fia_phynt_pa_z <- (fia_phynt_pa - mean(phynt_pa_null, na.rm=T))/sd(phynt_pa_null, na.rm=T)
		  
		  if (ncol(mat_p_noproblem) > 1) {
		  fia_funcpairwise_z <- (fia_funcpairwise - mean(funcpairwise_null, na.rm=T))/sd(funcpairwise_null, na.rm=T)
          fia_funcpairwise_pa_z <- (fia_funcpairwise_pa - mean(funcpairwise_pa_null, na.rm=T))/sd(funcpairwise_pa_null, na.rm=T)
          fia_funcnt_z <- (fia_funcnt - mean(funcnt_null, na.rm=T))/sd(funcnt_null, na.rm=T)
          fia_funcnt_pa_z <- (fia_funcnt_pa - mean(funcnt_pa_null, na.rm=T))/sd(funcnt_pa_null, na.rm=T)
		}
		else {
			fia_funcpairwise_z <- NA
			fia_funcpairwise_pa_z <- NA
			fia_funcnt_z <- NA
			fia_funcnt_pa_z <- NA
		}
			fia_list[[p]] <- data.frame(nneighb = nrow(mat_p) - 1, 
										beta_pd_pairwise_abundance = fia_phypairwise,
										beta_pd_pairwise_abundance_z = fia_phypairwise_z,
										beta_pd_nearest_abundance = fia_phynt,
										beta_pd_nearest_abundance_z = fia_phynt_z,
										beta_fd_pairwise_abundance = fia_funcpairwise,
										beta_fd_pairwise_abundance_z = fia_funcpairwise_z,
										beta_fd_nearest_abundance = fia_funcnt,
										beta_fd_nearest_abundance_z = fia_funcnt_z,
										beta_pd_pairwise_presence = fia_phypairwise_pa,
										beta_pd_pairwise_presence_z = fia_phypairwise_pa_z,
										beta_pd_nearest_presence = fia_phynt_pa,
										beta_pd_nearest_presence_z = fia_phynt_pa_z,
										beta_fd_pairwise_presence = fia_funcpairwise_pa,
										beta_fd_pairwise_presence_z = fia_funcpairwise_pa_z,
										beta_fd_nearest_presence = fia_funcnt_pa,
										beta_fd_nearest_presence_z = fia_funcnt_pa_z
			)
			
		  }
		  else {
			fia_list[[p]] <- data.frame(nneighb = NA, 
										beta_pd_pairwise_abundance = NA,
										beta_pd_pairwise_abundance_z = NA,
										beta_pd_nearest_abundance = NA,
										beta_pd_nearest_abundance_z = NA,
										beta_fd_pairwise_abundance = NA,
										beta_fd_pairwise_abundance_z = NA,
										beta_fd_nearest_abundance = NA,
										beta_fd_nearest_abundance_z = NA,
										beta_pd_pairwise_presence = NA,
										beta_pd_pairwise_presence_z = NA,
										beta_pd_nearest_presence = NA,
										beta_pd_nearest_presence_z = NA,
										beta_fd_pairwise_presence = NA,
										beta_fd_pairwise_presence_z = NA,
										beta_fd_nearest_presence = NA,
										beta_fd_nearest_presence_z = NA
			)
		  }
		}
		else {
		  fia_list[[p]] <- data.frame(nneighb = NA, 
										beta_pd_pairwise_abundance = NA,
										beta_pd_pairwise_abundance_z = NA,
										beta_pd_nearest_abundance = NA,
										beta_pd_nearest_abundance_z = NA,
										beta_fd_pairwise_abundance = NA,
										beta_fd_pairwise_abundance_z = NA,
										beta_fd_nearest_abundance = NA,
										beta_fd_nearest_abundance_z = NA,
										beta_pd_pairwise_presence = NA,
										beta_pd_pairwise_presence_z = NA,
										beta_pd_nearest_presence = NA,
										beta_pd_nearest_presence_z = NA,
										beta_fd_pairwise_presence = NA,
										beta_fd_pairwise_presence_z = NA,
										beta_fd_nearest_presence = NA,
										beta_fd_nearest_presence_z = NA
			)
		}
}

close(pb)

# Compile all of these values into a single data frame and save.
fia_betadiv <- do.call('rbind', fia_list)

write.csv(fia_betadiv, file = paste0('/mnt/research/nasabio/data/fia/betaoutput/tdpdfd_',r,'.csv'), row.names = FALSE)						  