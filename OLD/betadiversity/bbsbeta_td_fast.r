radii <- c(1000,2000,3000,4000,5000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
n_slices <- 50

combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]

load(paste0('/mnt/research/nasabio/data/bbs/mats/mat_', r, '_', slice, '.r'))
bbs_list <- list()

library(vegan)

pb <- txtProgressBar(0, length(all_mats), style=3)

for (p in 1:length(all_mats)) {
	setTxtProgressBar(pb, p)
	mat_p <- all_mats[[p]]
	
	if(inherits(mat_p, 'matrix')) {
		  if(nrow(mat_p) > 1 & ncol(mat_p) > 0) {
			
			dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
							
			bbs_list[[p]] <- data.frame(nneighb = nrow(mat_p) - 1, beta_td_pairwise_presence = mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard')))
			
		  }
		  else {
			bbs_list[[p]] <- data.frame(nneighb = NA, beta_td_pairwise_presence = NA)
		  }
		}
		else {
		  bbs_list[[p]] <- data.frame(nneighb = NA, beta_td_pairwise_presence = NA)
		}
}

close(pb)

# Compile all of these values into a single data frame and save.
bbs_betadiv <- do.call('rbind', bbs_list)

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/betaoutput/td_',r,'_',slice,'.csv'), row.names = FALSE)						  