load('C:/Users/Q/Dropbox/projects/nasabiodiv/fiaworkspace2.r')
source('betadiversity/pairwise_beta_focal.r')

library(sp)
library(vegan)
source('eda/fixpicante.r')


nnull <- 99
trydist <- as.matrix(trydist)

radii <- c(5, 10, 15, 20, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500)
n_slices <- 250
slice <- 51

gamma_div[p, r, ] <- diversity_3ways(m = neighbs, flavor = 'gamma', 
                                     td=T, pd=T, fd=T, abundance=T,
                                     pddist = fiadist, fddist = trydist,
                                     nnull = 5,
                                     phylo_spp = pnwphylo$tip.label, func_problem_spp = problemspp, phy = pnwphylo
)