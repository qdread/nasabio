# FIA
# Calculate alpha diversity for all plots.

# Split into 1500 groups.
# Use precalculated matrix.
# Loop through each FIA plot.
# Calculate all alpha diversity.

# Edited 21 Dec 2017 to use the whole USA

load('/mnt/research/nasabio/data/fia/fiaworkspace_nospatial_wholeusa.r')
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa.csv', stringsAsFactors = FALSE)
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
source('/mnt/research/nasabio/code/nofuncspp.r')

library(sp)
library(vegan)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 99
trydist <- as.matrix(trydist)

n_slices <- 1500
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaplotmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

alpha_div <- diversity_3ways(m = fiaplotmat[rowidxmin:rowidxmax,], flavor = 'alpha', 
											 dotd = TRUE, dopd = TRUE, dofd = TRUE, abundance = TRUE,
											 pddist = fiadist, fddist = trydist,
											 nnull = nnull,
											 phylo_spp = fullphylo$tip.label, func_problem_spp = nofuncspp, combine = FALSE)

alpha_div <- cbind(plotmetadata[rowidxmin:rowidxmax,], alpha_div)
save(alpha_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/usa/alpha_', slice, '.r'))
