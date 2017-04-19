radii <- c(50000, 75000, 100000, 150000, 200000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

r <- radii[task]


load(paste0('/mnt/research/nasabio/data/bbs/mats/routemat_', r, '.r'))
bbs_list <- list()

load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(vegan)
library(vegetarian)
source('~/code/fia/fixpicante.r')

nnull <- 99 # Reduce to save time

pb <- txtProgressBar(0, length(all_mats), style=3)

for (p in 1:length(all_mats)) {
	setTxtProgressBar(pb, p)
	mat_p <- all_mats[[p]]
	
	if(inherits(mat_p, 'matrix')) {
		  if(nrow(mat_p) > 1 & ncol(mat_p) > 1) {
			
			dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
			
			# Observed statistics
			beta_td_pairwise_presence = mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
			bbs_phypairwise_pa <- mean(comdist(comm = mat_p, dis = ericdist, abundance.weighted = FALSE))
            bbs_phynt_pa <- mean(comdistnt(comm = mat_p, dis = ericdist, abundance.weighted = FALSE))
            bbs_funcpairwise_pa <- mean(comdist(comm = mat_p, dis = birdtraitdist, abundance.weighted = FALSE))
            bbs_funcnt_pa <- mean(comdistnt(comm = mat_p, dis = birdtraitdist, abundance.weighted = FALSE))
			
			# Null distribution of statistics
			phypairwise_pa_null <- phynt_pa_null <- rep(NA, nnull)
		    funcpairwise_pa_null <- funcnt_pa_null <- rep(NA, nnull)
          
          for (sim in 1:nnull) {
            ericnullidx <- sample(1:nrow(ericdist))
            ericdistnull <- ericdist
            dimnames(ericdistnull)[[1]] <- dimnames(ericdistnull)[[1]][ericnullidx]
            dimnames(ericdistnull)[[2]] <- dimnames(ericdistnull)[[2]][ericnullidx]
            
			birdtraitnullidx <- sample(1:nrow(birdtraitdist))
            birdtraitdistnull <- birdtraitdist
            dimnames(birdtraitdistnull)[[1]] <- dimnames(birdtraitdistnull)[[1]][birdtraitnullidx]
            dimnames(birdtraitdistnull)[[2]] <- dimnames(birdtraitdistnull)[[2]][birdtraitnullidx]
			
            phypairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = ericdistnull, abundance.weighted = FALSE))
            phynt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = ericdistnull, abundance.weighted = FALSE))
			
            funcpairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = birdtraitdistnull, abundance.weighted = FALSE))
            funcnt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = birdtraitdistnull, abundance.weighted = FALSE))
		  }
		  
		  bbs_phypairwise_pa_z <- (bbs_phypairwise_pa - mean(phypairwise_pa_null, na.rm=T))/sd(phypairwise_pa_null, na.rm=T)
          bbs_phynt_pa_z <- (bbs_phynt_pa - mean(phynt_pa_null, na.rm=T))/sd(phynt_pa_null, na.rm=T)
		  
          bbs_funcpairwise_pa_z <- (bbs_funcpairwise_pa - mean(funcpairwise_pa_null, na.rm=T))/sd(funcpairwise_pa_null, na.rm=T)
          bbs_funcnt_pa_z <- (bbs_funcnt_pa - mean(funcnt_pa_null, na.rm=T))/sd(funcnt_pa_null, na.rm=T)

			
			bbs_list[[p]] <- data.frame(nneighb = nrow(mat_p) - 1, 
										beta_td_pairwise_presence = beta_td_pairwise_presence,
										beta_pd_pairwise_presence = bbs_phypairwise_pa,
										beta_pd_pairwise_presence_z = bbs_phypairwise_pa_z,
										beta_pd_nearest_presence = bbs_phynt_pa,
										beta_pd_nearest_presence_z = bbs_phynt_pa_z,
										beta_fd_pairwise_presence = bbs_funcpairwise_pa,
										beta_fd_pairwise_presence_z = bbs_funcpairwise_pa_z,
										beta_fd_nearest_presence = bbs_funcnt_pa,
										beta_fd_nearest_presence_z = bbs_funcnt_pa_z
			)
			
		  }
		  else {
			bbs_list[[p]] <- data.frame(nneighb = NA, 
										beta_td_pairwise_presence = NA,
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
		  bbs_list[[p]] <- data.frame(nneighb = NA, 
										beta_td_pairwise_presence = NA,
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
bbs_betadiv <- do.call('rbind', bbs_list)

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/betaoutput/route_tdpdfd_',r,'.csv'), row.names = FALSE)						  