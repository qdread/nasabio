# Created on 09 June 2017 to run beta diversity using Baselga and Podani methods (taxonomic only)
# Edit 16 Aug: very small slices
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
radii <- c(1, 5, 7.5, 10, 20, 50, 75, 100, 150, 200, 300) * 1000
n_slices <- 250


combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]


load(paste0('/mnt/research/nasabio/data/fia/mats/newmat_', as.character(as.integer(r)), '.r'))
fia_list <- list()

load('/mnt/research/nasabio/data/fia/fiaworkspace2.r') 

source('~/code/fia/beta_part_finalindex.r')

null_result <- data.frame(index = rep(c('sorensen','bray','sorensen'), each=3),
                          diversity = rep(c('taxonomic','phylogenetic'), times=c(6,3)),
                          partition = c('replacement', 'nestedness', 'total'),
                          abundance = rep(c(FALSE, TRUE, FALSE), each=3),
                          beta = NA)

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
        beta_p <- beta_part(m = mat_p, abundance = TRUE, pairwise = FALSE, index_family = 'sorensen', TD=TRUE, PD=TRUE, FD=FALSE, phylo_tree = pnwphylo)
	   
			fia_list[[p]] <- cbind(nneighb = nrow(mat_p) - 1, beta_p)
	
		  }
		  else {
			fia_list[[p]] <- cbind(nneighb = NA, null_result)
		  }
		}
		else {
		  fia_list[[p]] <- cbind(nneighb = NA, null_result)
		}
		}

close(pb)

# Compile all of these values into a single data frame and save.
fia_betadiv <- as.data.frame(do.call('rbind', fia_list))

write.csv(fia_betadiv, file = paste0('/mnt/research/nasabio/data/fia/betaoutput/sliced/betapartfinal_',as.character(as.integer(r)),'_',slice,'.csv'), row.names = FALSE)						  
