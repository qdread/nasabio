task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
radii <- c(50, 75, 100, 150, 200, 300) * 1000
n_slices <- 20

# Created on 9 June to run new beta diversity.

combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]


load(paste0('/mnt/research/nasabio/data/bbs/mats/newroutemat_', as.character(as.integer(r)), '.r'))
bbs_list <- list()

# Determine row indices for the slice of the matrix to be used.
idx <- round(seq(0,length(all_mats),length.out=n_slices + 1))
idxmin <- idx[slice]+1
idxmax <- idx[slice+1]

source('~/code/fia/beta_part.r')

null_result <- 	res <- data.frame(family = rep(c('podani', 'baselga'), each = 10),
					  index = rep(c('sorensen', 'jaccard'), each = 5),
					  divtype = c('total','replacement','nestedness','replacement_proportion','nestedness_proportion'),
					  abundance = FALSE,
					  beta = NA)

pb <- txtProgressBar(0, length(idxmin:idxmax), style=3)

for (p in idxmin:idxmax) {
	setTxtProgressBar(pb, p)
	mat_p <- all_mats[[p]]
	
	if (inherits(mat_p, 'matrix')) {
		  if (nrow(mat_p) > 1 & ncol(mat_p) > 1) {
			
			dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
			
			beta_td <- beta_baselga_podani(m = mat_p, abundance = FALSE)
			
			bbs_list[[(p - idxmin + 1)]] <- cbind(nneighb = nrow(mat_p) - 1, beta_td)
			
		  }
		  else {
			bbs_list[[(p - idxmin + 1)]] <- cbind(nneighb = NA, null_result)
		  }
		}
		else {
		  bbs_list[[(p - idxmin + 1)]] <- cbind(nneighb = NA, null_result)
		}
}

close(pb)

# Compile all of these values into a single data frame and save.
bbs_betadiv <- do.call('rbind', bbs_list)

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/betaoutput/byroute/route_betabaselga_',as.character(as.integer(r)),'_',slice,'.csv'), row.names = FALSE)						  
