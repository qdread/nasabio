# Calculate gamma diversity within different radii for BBS.
# in contrast to FIA, we can only use groups within the same year. 1997-present.

# Alpha diversity for all together since it is not necessary to do parallel (update 31 Jan)

load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

# Replace AOU codes in the trait matrix with actual species names.

library(sp)
library(vegan)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 999 # Edited 12 Dec 2018

alpha_div <- diversity_3ways(m = bbsmat_byroute_oneyear, flavor = 'alpha', 
                             dotd=T, dopd=T, dofd=T, abundance=F,
                             pddist = ericdist, fddist = birdtraitdist,
                             nnull = nnull,
                             phylo_spp = NULL, func_problem_spp = NULL, combine = FALSE)

bbs_alphadiv <- cbind(bbscov_oneyear, alpha_div)
write.csv(bbs_alphadiv, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/bbs_alphadiv.csv', row.names = FALSE)
