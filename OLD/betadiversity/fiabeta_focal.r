# fork created 15 may 2017

radii <- c(10000,20000,50000,100000)
n_slices <- 50
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]


load(paste0('/mnt/research/nasabio/data/fia/mats/newmat_', as.character(as.integer(r)), '.r'))
fia_list <- list()

load('/mnt/research/nasabio/data/fia/fiaworkspace.r') 

library(vegan)
library(vegetarian)
source('~/code/fia/fixpicante.r')
source('~/code/fia/pairwise_beta_focal.r')

nnull <- 99 # Reduce to save time

trydist <- as.matrix(trydist)

# Determine row indices for the slice of the matrix to be used.
idx <- round(seq(0,length(all_mats),length.out=n_slices + 1))
idxmin <- idx[slice]+1
idxmax <- idx[slice+1]

pb <- txtProgressBar(0, length(all_mats), style=3)

for (p in 1:length(idxmin:idxmax)) {
	setTxtProgressBar(pb, p)
	mat_p <- all_mats[[(p + idxmin - 1)]]
	
	if(inherits(mat_p, 'matrix')) {
		  if(nrow(mat_p) > 1 & ncol(mat_p) > 1) {
			
			# Fix the species names to match the phylogeny, and get rid of the unknown species.
        
        mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% pnwphylo$tip.label, drop = FALSE]
        mat_p_noproblem <- mat_p[, !dimnames(mat_p)[[2]] %in% problemspp, drop = FALSE]
        
        # Calculate beta-diversity for that matrix.
        fia_betatdpd <- pairwise_beta_focal(m=mat_p, td=T, pd=T, fd=F, abundance=T, pddist=fiadist)

		  if (ncol(mat_p_noproblem) > 1) {
			  fia_betafd <- pairwise_beta_focal(m=mat_p, td=F, pd=F, fd=T, abundance=T, fddist=trydist)
		  }
		  else {
			  fia_betafd <- c(beta_fd_pairwise_pa=NA, beta_fd_nt_pa=NA, beta_fd_pairwise=NA, beta_fd_nt=NA)	
		  }
          # Null models by scrambling distance matrix
		  
          phy_null <- matrix(NA, nrow = nnull, ncol = 4)
		  func_null <- matrix(NA, nrow = nnull, ncol = 4)
          
          for (sim in 1:nnull) {
            nullidx <- sample(1:nrow(fiadist))
            fiadistnull <- fiadist
            dimnames(fiadistnull)[[1]] <- dimnames(fiadistnull)[[1]][nullidx]
            dimnames(fiadistnull)[[2]] <- dimnames(fiadistnull)[[2]][nullidx]
            
			trynullidx <- sample(1:nrow(trydist))
            trydistnull <- trydist
            dimnames(trydistnull)[[1]] <- dimnames(trydistnull)[[1]][trynullidx]
            dimnames(trydistnull)[[2]] <- dimnames(trydistnull)[[2]][trynullidx]
			
			phy_null[sim, ] <- pairwise_beta_focal(m=mat_p, td=F, pd=T, fd=F, abundance=T, pddist=fiadistnull)
			
			if (ncol(mat_p_noproblem) > 1) {
				func_null[sim, ] <- pairwise_beta_focal(m=mat_p, td=F, pd=F, fd=T, abundance=T, fddist=trydistnull)
			}
			
          }
          
		  phypairwise_pa_z <- (fia_betatdpd['beta_pd_pairwise_pa'] - mean(phy_null[,1], na.rm=T))/sd(phy_null[,1], na.rm=T)
          phynt_pa_z <- (fia_betatdpd['beta_pd_nt_pa'] - mean(phy_null[,2], na.rm=T))/sd(phy_null[,2], na.rm=T)
		  phypairwise_z <- (fia_betatdpd['beta_pd_pairwise'] - mean(phy_null[,3], na.rm=T))/sd(phy_null[,3], na.rm=T)
          phynt_z <- (fia_betatdpd['beta_pd_nt'] - mean(phy_null[,4], na.rm=T))/sd(phy_null[,4], na.rm=T)
		  
		  if (ncol(mat_p_noproblem) > 1) {
			funcpairwise_pa_z <- (fia_betafd['beta_fd_pairwise_pa'] - mean(func_null[,1], na.rm=T))/sd(func_null[,1], na.rm=T)
            funcnt_pa_z <- (fia_betafd['beta_fd_nt_pa'] - mean(func_null[,2], na.rm=T))/sd(func_null[,2], na.rm=T)
			funcpairwise_pa_z <- (fia_betafd['beta_fd_pairwise'] - mean(func_null[,3], na.rm=T))/sd(func_null[,3], na.rm=T)
            funcnt_pa_z <- (fia_betafd['beta_fd_nt'] - mean(func_null[,4], na.rm=T))/sd(func_null[,4], na.rm=T)
		  }
		else {
			func_null_z <- rep(NA, 4)
		}
			fia_list[[p]] <- c(nneighb = nrow(mat_p) - 1, 
										fia_betatdpd,
										fia_betafd,
										phypairwise_pa_z, phynt_pa_z, phypairwise_z, phynt_z,
										funcpairwise_pa_z, funcnt_pa_z, funcpairwise_pa_z, funcnt_pa_z)
	
		  }
		  else {
			fia_list[[p]] <- rep(NA,20)
		  }
		}
		else {
		  fia_list[[p]] <- rep(NA,20)
		}
		names(fia_list[[p]]) <- c('nneighb', 'beta_td_pairwise_presence', 'beta_td_shannon', 'beta_td_pairwise',
								  'beta_pd_pairwise_presence', 'beta_pd_nt_presence', 'beta_pd_pairwise', 'beta_pd_nt',
								  'beta_fd_pairwise_presence', 'beta_fd_nt_presence', 'beta_fd_pairwise', 'beta_fd_nt',
								  'beta_pd_pairwise_presence_z', 'beta_pd_nt_presence_z', 'beta_pd_pairwise_z', 'beta_pd_nt_z',
								  'beta_fd_pairwise_presence_z', 'beta_fd_nt_presence_z', 'beta_fd_pairwise_z', 'beta_fd_nt_z')
}

close(pb)

# Compile all of these values into a single data frame and save.
fia_betadiv <- as.data.frame(do.call('rbind', fia_list))

write.csv(fia_betadiv, file = paste0('/mnt/research/nasabio/data/fia/betaoutput/sliced/focal_tdpdfd_',as.character(as.integer(r)),'_',slice,'.csv'), row.names = FALSE)						  