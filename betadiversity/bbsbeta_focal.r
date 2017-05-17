# fork created 15 may 2017

radii <- c(50000, 75000, 100000, 150000, 200000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
n_slices <- 50

combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]


load(paste0('/mnt/research/nasabio/data/bbs/mats/newroutemat_', as.character(as.integer(r)), '.r'))
bbs_list <- list()

# Determine row indices for the slice of the matrix to be used.
idx <- round(seq(0,length(all_mats),length.out=n_slices + 1))
idxmin <- idx[slice]+1
idxmax <- idx[slice+1]

load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(vegan)
library(vegetarian)
source('~/code/fia/fixpicante.r')
source('~/code/fia/pairwise_beta_focal.r')

nnull <- 99 # Reduce to save time

pb <- txtProgressBar(0, length(idxmin:idxmax), style=3)

for (p in 1:length(idxmin:idxmax)) {
	setTxtProgressBar(pb, p)
	mat_p <- all_mats[[(p + idxmin - 1)]]
	
	if(inherits(mat_p, 'matrix')) {
		  if(nrow(mat_p) > 1 & ncol(mat_p) > 1) {
			
			dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
			
			# Observed statistics
			bbs_betatdpdfd <- pairwise_beta_focal(m=mat_p, td=T, pd=T, fd=T, abundance=F, pddist=ericdist, fddist=birdtraitdist)
			
			# Null distribution of statistics
			phy_null <- matrix(NA, nrow = nnull, ncol = 2)
		    func_null <- matrix(NA, nrow = nnull, ncol = 2)
          
          for (sim in 1:nnull) {
            ericnullidx <- sample(1:nrow(ericdist))
            ericdistnull <- ericdist
            dimnames(ericdistnull)[[1]] <- dimnames(ericdistnull)[[1]][ericnullidx]
            dimnames(ericdistnull)[[2]] <- dimnames(ericdistnull)[[2]][ericnullidx]
            
			birdtraitnullidx <- sample(1:nrow(birdtraitdist))
            birdtraitdistnull <- birdtraitdist
            dimnames(birdtraitdistnull)[[1]] <- dimnames(birdtraitdistnull)[[1]][birdtraitnullidx]
            dimnames(birdtraitdistnull)[[2]] <- dimnames(birdtraitdistnull)[[2]][birdtraitnullidx]
			
			phy_null[sim, ] <- pairwise_beta_focal(m=mat_p, td=F, pd=T, fd=F, abundance=F, pddist=ericdistnull)
			func_null[sim, ] <- pairwise_beta_focal(m=mat_p, td=F, pd=F, fd=T, abundance=F, fddist=birdtraitdistnull)
		  }
		  
		  phypairwise_pa_z <- (bbs_betatdpdfd['beta_pd_pairwise_pa'] - mean(phy_null[,1], na.rm=T))/sd(phy_null[,1], na.rm=T)
          phynt_pa_z <- (bbs_betatdpdfd['beta_pd_nt_pa'] - mean(phy_null[,2], na.rm=T))/sd(phy_null[,2], na.rm=T)
		  
          funcpairwise_pa_z <- (bbs_betatdpdfd['beta_fd_pairwise_pa'] - mean(func_null[,1], na.rm=T))/sd(func_null[,1], na.rm=T)
          funcnt_pa_z <- (bbs_betatdpdfd['beta_fd_nt_pa'] - mean(func_null[,2], na.rm=T))/sd(func_null[,2], na.rm=T)
			
			bbs_list[[p]] <- c(nneighb = nrow(mat_p) - 1, 
										bbs_betatdpdfd,
										phypairwise_pa_z, phynt_pa_z,
										funcpairwise_pa_z, funcnt_pa_z)
			
		  }
		  else {
			bbs_list[[p]] <- rep(NA,10)
		  }
		}
		else {
		  bbs_list[[p]] <- rep(NA,10)
		}
		names(bbs_list[[p]]) <- c('nneighb', 'beta_td_pairwise_presence', 
								  'beta_pd_pairwise_presence', 'beta_pd_nt_presence', 
								  'beta_fd_pairwise_presence', 'beta_fd_nt_presence', 
								  'beta_pd_pairwise_presence_z', 'beta_pd_nt_presence_z',
								  'beta_fd_pairwise_presence_z', 'beta_fd_nt_presence_z')
}

close(pb)

# Compile all of these values into a single data frame and save.
bbs_betadiv <- as.data.frame(do.call('rbind', bbs_list))

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/betaoutput/byroute/focal_tdpdfd_',as.character(as.integer(r)),'_',slice,'.csv'), row.names = FALSE)						  
