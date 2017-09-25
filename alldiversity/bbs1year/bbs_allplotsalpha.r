# Calculate gamma diversity within different radii for BBS.
# in contrast to FIA, we can only use groups within the same year. 1997-present.

# Split into 250 groups.
# Use precalculated matrix.
# Loop through each BBS plot.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 500 km of the target plot.
# Do median alpha, and total gamma, for all the neighbors within each radius.
# All others outside that radius get NA.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')
source('~/code/fia/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

# Replace AOU codes in the trait matrix with actual species names.

library(sp)
library(vegan)
source('~/code/fia/fixpicante.r')

nnull <- 99

n_slices <- 250
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(bbsmat_byroute_oneyear),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

alpha_div <- diversity_3ways(m = bbsmat_byroute_oneyear[rowidxmin:rowidxmax,], flavor = 'alpha', 
                             dotd=T, dopd=T, dofd=T, abundance=F,
                             pddist = ericdist, fddist = birdtraitdist,
                             nnull = nnull,
                             phylo_spp = NULL, func_problem_spp = NULL, combine = F)

alpha_div <- cbind(bbscov_oneyear[rowidxmin:rowidxmax,], alpha_div)
save(alpha_div, file = paste0('/mnt/research/nasabio/data/bbs/diversity1year/alpha_', slice, '.r'))
