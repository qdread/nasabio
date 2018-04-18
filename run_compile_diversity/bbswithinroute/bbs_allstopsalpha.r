# Calculate gamma diversity within different radii for BBS.
# in contrast to FIA, we can only use groups within the same year. 1997-present.

# Alpha diversity for all together since it is not necessary to do parallel
# In another script, load this output and take the averages by the different route subsets.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_bystop_20072016.r')
source('~/code/fia/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(sp)
library(vegan)

nnull <- 99

# Get the alpha diversity for each stop, then combine to 

alpha_div <- diversity_3ways(m = bbsmat_oneyear, flavor = 'alpha', 
                             dotd=T, dopd=T, dofd=T, abundance=F,
                             pddist = ericdist, fddist = birdtraitdist,
                             nnull = nnull,
                             phylo_spp = NULL, func_problem_spp = NULL, combine = FALSE)

bbs_alphadiv <- cbind(bbscov_oneyear, alpha_div)
write.csv(bbs_alphadiv, file = '/mnt/research/nasabio/data/bbs/diversitywithinroute/bbs_alphadiv.csv', row.names = FALSE)
