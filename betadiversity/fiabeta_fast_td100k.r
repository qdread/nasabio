n_slices <- 50
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
r <- 1e5

load(paste0('/mnt/research/nasabio/data/fia/mats/mat_', as.character(as.integer(r)), '.r'))
fia_list <- list()

load('/mnt/research/nasabio/data/fia/fiaworkspace.r') 

 library(vegan)
 library(vegetarian)

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
       # mat_p_noproblem <- mat_p[, !dimnames(mat_p)[[2]] %in% problemspp, drop = FALSE]
        
        # Calculate beta-diversity for that matrix.
        
        fia_shannonbetadiv <- d(abundances = mat_p, lev = 'beta', wts = FALSE, q = 1)
        fia_meanpairwisedissim <- mean(vegdist(x = mat_p, binary = FALSE, method = 'jaccard'), na.rm=TRUE)
        fia_meanpairwisedissim_pa <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'), na.rm=TRUE)
        
			fia_list[[p]] <- data.frame(nneighb = nrow(mat_p) - 1, 
										beta_shannon = fia_shannonbetadiv,
										beta_pairwise_abundance = fia_meanpairwisedissim,
										beta_pairwise_presence = fia_meanpairwisedissim_pa
										
			)
			
		  }
		  else {
			fia_list[[p]] <- data.frame(nneighb = NA,
										beta_shannon = NA,
										beta_pairwise_abundance = NA,
										beta_pairwise_presence = NA
										
			)
		  }
		}
		else {
		  fia_list[[p]] <- data.frame(nneighb = NA,
										beta_shannon = NA,
										beta_pairwise_abundance = NA,
										beta_pairwise_presence = NA
										
			)
		}
}

close(pb)

# Compile all of these values into a single data frame and save.
fia_betadiv <- do.call('rbind', fia_list)

write.csv(fia_betadiv, file = paste0('/mnt/research/nasabio/data/fia/betaoutput/sliced/td_',as.character(as.integer(r)),'_',slice,'.csv'), row.names = FALSE)						  