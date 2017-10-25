# in each slice, load calculated fia matrix.
# Calculate gamma diversity within different radii.

# Split into 250 groups.
# Use precalculated matrix.
# Loop through each FIA plot.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 500 km of the target plot.
# Do median alpha, and total gamma, for all the neighbors within each radius.
# All others outside that radius get NA.

load('/mnt/research/nasabio/data/fia/fiaworkspace_nospatial.r')
source('~/code/fia/pairwise_beta_focal.r')

library(sp)
library(vegan)
source('~/code/fia/fixpicante.r')

nnull <- 99
trydist <- as.matrix(trydist)

radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300)
n_slices <- 250
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaplotmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Declare structures to hold data

pb <- txtProgressBar(0, length(radii), style = 3)
gamma_div <- array(NA, dim = c(nrow(fiaplotmat), length(radii), 11))

cnames <- c('richness', 'shannon', 'evenness',
            'MPD_pa_z', 'MNTD_pa_z',
            'MPD_z', 'MNTD_z',
            'MPDfunc_pa_z', 'MNTDfunc_pa_z',
            'MPDfunc_z', 'MNTDfunc_z')

null_result <- rep(NA, length(cnames))

for (r in 1:length(radii)) {
	load(paste0('/mnt/research/nasabio/data/fia/mats/unfuzzedmat_', as.character(as.integer(radii[r] * 1000)), '.r'))
	setTxtProgressBar(pb, r)
	for (p in rowidxmin:rowidxmax) {
		neighbs <- all_mats[[p]]
		gamma_div[p, r, ] <- tryCatch(diversity_3ways(m = neighbs, flavor = 'gamma', 
											 dotd=T, dopd=T, dofd=T, abundance=T,
											 pddist = fiadist, fddist = trydist,
											 nnull = nnull,
											 phylo_spp = pnwphylo$tip.label, func_problem_spp = NULL),
									  error = function(e) null_result)	
	}

}

close(pb)



dimnames(gamma_div)[[3]] <- cnames
dimnames(gamma_div)[[2]] <- paste('r',radii,sep='_')

save(gamma_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/unfuzzed/gamma_', slice, '.r'))
