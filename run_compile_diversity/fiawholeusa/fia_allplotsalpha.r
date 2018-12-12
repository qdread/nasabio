# FIA
# Calculate alpha diversity for all plots.

# Split into 1000 groups.
# Use precalculated matrix.
# Loop through each FIA plot.
# Calculate all alpha diversity.

# Edited 21 Dec 2017 to use the whole USA
# Edited 26 Nov 2018: file paths update, update for slurm
# Edited 27 Nov 2018: get rid of parsing
# edited 12 Dec 2018: increase null model iterations

load('/mnt/ffs17/groups/nasabio/fiaworkspace_nospatial_wholeusa_2018.r')
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa_2018.csv', stringsAsFactors = FALSE)
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
source('/mnt/research/nasabio/code/nofuncspp.r')

library(sp)
library(vegan)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 999
trydist <- as.matrix(trydist)

n_slices <- 1000
slice <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaplotmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

alpha_div <- diversity_3ways(m = fiaplotmat[rowidxmin:rowidxmax,], flavor = 'alpha', 
											 dotd = TRUE, dopd = TRUE, dofd = TRUE, abundance = TRUE,
											 pddist = fiadist, fddist = trydist,
											 nnull = nnull,
											 phylo_spp = fullphylo$tip.label, func_problem_spp = nofuncspp, combine = FALSE)

alpha_div <- cbind(PLT_CN = plotmetadata[rowidxmin:rowidxmax,], alpha_div)
save(alpha_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/alpha_', slice, '.r'))
