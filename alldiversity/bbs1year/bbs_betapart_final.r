# Redo BBS beta-partition with "final" index.
# This will include taxonomic and phylogenetic diversity.

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
radii <- c(50, 75, 100, 150, 200, 300, 400, 500) * 1000
n_slices <- 250

combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]

# Load matrix
load(paste0('/mnt/research/nasabio/data/bbs/mats/oneyearroutemat_', as.character(as.integer(r)), '.r'))

# Must load bbs phylogenetic tree here.
library(ape)
erictree <- read.tree('/mnt/research/nasabio/data/bbs/ericson_cons.tre')

# Load species list and "fix" the consensus tree to have correct labels that match the labels on the community matrix.
bbsspp <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/specieslist.csv', stringsAsFactors = FALSE)

tlabel1 <- erictree$tip.label
tlabel1 <- gsub('_', ' ', tlabel1)

phymatchidx <- rep(NA,length(tlabel1))

for (i in 1:length(tlabel1)) {
	phymatchidx[i] <- c(which(bbsspp$Latin_Name_clean == tlabel1[i] | bbsspp$Latin_Name_synonym == tlabel1[i] | bbsspp$Latin_Name_synonym2 == tlabel1[i]), NA)[1]
}
tlabelaou <- bbsspp$AOU[phymatchidx]

dimnames_tlabel <- rep(NA, length(tlabelaou))
for (i in 1:length(tlabelaou)) {
	names_i <- bbsspp[bbsspp$AOU == tlabelaou[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_tlabel[i] <- names_i[1]
}

erictree$tip.label <- dimnames_tlabel

# Source the beta partitioning functions.
source('~/code/fia/beta_part_finalindex.r')

null_result <- data.frame(index = 'sorensen',
                          diversity = rep(c('taxonomic','phylogenetic'), each = 3),
                          partition = c('replacement', 'nestedness', 'total'),
                          abundance = FALSE,
                          beta = NA)

# Determine row indices for the slice of the matrix to be used.
idx <- round(seq(0,length(all_mats),length.out=n_slices + 1))
idxmin <- idx[slice]+1
idxmax <- idx[slice+1]

bbs_list <- list()

pb <- txtProgressBar(0, length(all_mats), style=3)

for (p in 1:length(idxmin:idxmax)) {
	setTxtProgressBar(pb, p)
	mat_p <- all_mats[[(p + idxmin - 1)]]
	
	if(inherits(mat_p, 'matrix')) {
		  if(nrow(mat_p) > 1 & ncol(mat_p) > 1) {
			
			# Fix the species names to match the phylogeny, and get rid of the unknown species.
        beta_p <- beta_part(m = mat_p, abundance = FALSE, pairwise = FALSE, index_family = 'sorensen', TD=TRUE, PD=TRUE, FD=FALSE, phylo_tree = erictree)
	   
			bbs_list[[p]] <- cbind(nneighb = nrow(mat_p) - 1, beta_p)
	
		  }
		  else {
			bbs_list[[p]] <- cbind(nneighb = NA, null_result)
		  }
	}
	else {
		  bbs_list[[p]] <- cbind(nneighb = NA, null_result)
	}
}

close(pb)

# Compile all of these values into a single data frame and save.
bbs_betadiv <- as.data.frame(do.call('rbind', bbs_list))

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/diversity1year/betapartfinal_',as.character(as.integer(r)),'_',slice,'.csv'), row.names = FALSE)				