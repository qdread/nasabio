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

p <- 9002
r <- 4

dist_p <- spDistsN1(pts=with(fiacoords, cbind(lon, lat)), pt = c(fiacoords$lon[p], fiacoords$lat[p]), longlat = TRUE)
neighbs <- fiaplotmat[dist_p <= radii[r], ]
gamma_div[p, r, ] <- diversity_3ways(m = neighbs, flavor = 'gamma', 
                                     dotd=T, dopd=T, dofd=T, abundance=T,
                                     pddist = fiadist, fddist = trydist,
                                     nnull = 99,
                                     phylo_spp = pnwphylo$tip.label, func_problem_spp = problemspp, phy = pnwphylo
)