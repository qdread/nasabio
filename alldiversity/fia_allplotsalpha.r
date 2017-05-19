# FIA
# Calculate alpha diversity for all plots.

# Split into 250 groups.
# Use precalculated matrix.
# Loop through each FIA plot.
# Calculate all alpha diversity.

load('/mnt/research/nasabio/data/fia/fiaworkspace2.r')
source('~/code/fia/pairwise_beta_focal.r')

library(sp)
library(vegan)
source('~/code/fia/fixpicante.r')

nnull <- 99
trydist <- as.matrix(trydist)

n_slices <- 250
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaplotmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

alpha_div <- diversity_3ways(m = fiaplotmat[rowidxmin:rowidxmax,], flavor = 'alpha', 
											 dotd=T, dopd=T, dofd=T, abundance=T,
											 pddist = fiadist, fddist = trydist,
											 nnull = nnull,
											 phylo_spp = pnwphylo$tip.label, func_problem_spp = problemspp, combine = F)

alpha_div <- cbind(fiacoords[rowidxmin:rowidxmax,], alpha_div)
save(alpha_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/alpha_', slice, '.r'))