load('C:/Users/Q/Dropbox/projects/nasabiodiv/fiaworkspace2.r')
source('betadiversity/pairwise_beta_focal.r')
load('C:/Users/Q/Dropbox/projects/nasabiodiv/trymat_clean.r')

library(sp)
library(vegan)
source('eda/fixpicante.r')


nnull <- 99
trydist <- as.matrix(trydist)

radii <- c(5, 10, 15, 20, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500)
n_slices <- 250
slice <- 51

p <- 4507
r <- 1

dist_p <- spDistsN1(pts=with(fiacoords, cbind(lon, lat)), pt = c(fiacoords$lon[p], fiacoords$lat[p]), longlat = TRUE)
neighbs <- fiaplotmat[dist_p <= radii[r], , drop = FALSE]

m <- neighbs
flavor <- 'gamma'
dotd <- T
dopd <- T
dofd <- T
abundance <- T
pddist <- fiadist
fddist <- trydist
fdmat <- try_noproblem
phylo_spp <- pnwphylo$tip.label
func_problem_spp <- problemspp

diversity_3ways(m = neighbs, flavor = 'gamma', 
                                     dotd=T, dopd=T, dofd=F, abundance=T,
                                     pddist = fiadist, fdmat = try_noproblem,
                                     nnull = 10,
                                     phylo_spp = pnwphylo$tip.label, func_problem_spp = problemspp)

############################################3
# bbs

load('C:/Users/Q/Dropbox/projects/nasabiodiv/bbsworkspace_byroute.r')
source('betadiversity/pairwise_beta_focal.r')
load('C:/Users/Q/Dropbox/projects/nasabiodiv/bbspdfddist.r') # Phy and Func distance matrices.
load('C:/Users/Q/Dropbox/projects/nasabiodiv/birdtraitmat_clean.r')

library(sp)
library(vegan)
source('eda/fixpicante.r')

nnull <- 10

# run through all radii in each task, but split by year and also slice it up some.
radii <- c(50, 75, 100, 150, 200, 250, 300, 400, 500)
years <- 1997:2015 # 19 years, so we can have 13 slices per year. 19*13 = 247
n_slices <- 13

year <- 2001
slice <- 10


# Split the master bbs matrix by year.
bbs_year_mat <- fixedbbsmat_byroute[bbscov$year == year, ]
bbs_year_coords <- bbscov[bbscov$year == year, ]

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(bbs_year_mat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

p <- 1580
r <- 3

dist_p <- spDistsN1(pts=with(bbs_year_coords, cbind(lon, lat)), pt = c(bbs_year_coords$lon[p], bbs_year_coords$lat[p]), longlat = TRUE)
neighbs <- bbs_year_mat[dist_p <= radii[r], , drop = FALSE]

m <- neighbs
flavor <- 'gamma'
dotd <- T
dopd <- T
dofd <- T
abundance <- F
pddist <- ericdist
fddist <- birdtraitdist
fdmat <- birdtraitclean
phylo_spp <- NULL
func_problem_spp <- NULL

diversity_3ways(m = neighbs, flavor = 'gamma', 
                dotd=T, dopd=T, dofd=T, abundance=F,
                pddist = ericdist, fdmat = birdtraitclean,
                nnull = nnull,
                phylo_spp = NULL, func_problem_spp = NULL)

