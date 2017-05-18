# in each slice, load calculated fia matrix.
# Calculate gamma diversity within different radii.

# Split into 250 groups.
# Use precalculated matrix.
# Loop through each FIA plot.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 500 km of the target plot.
# Do median alpha, and total gamma, for all the neighbors within each radius.
# All others outside that radius get NA.

load('/mnt/research/nasabio/data/fia/fiaworkspace2.r')
source('~/code/fia/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/fia/trymat_clean.r') # trait matrix.

library(sp)
library(vegan)
source('~/code/fia/fixpicante.r')

nnull <- 99
trydist <- as.matrix(trydist)

radii <- c(5, 10, 15, 20, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500)
n_slices <- 250
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaplotmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Declare structures to hold data

pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
gamma_div <- array(NA, dim = c(nrow(fiaplotmat), length(radii), 15))

for (p in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, p)
	# Distances between target plot and all other plots.
	dist_p <- spDistsN1(pts=with(fiacoords, cbind(lon, lat)), pt = c(fiacoords$lon[p], fiacoords$lat[p]), longlat = TRUE)
	
	for (r in 1:length(radii)) {
		neighbs <- fiaplotmat[dist_p <= radii[r], , drop = FALSE]
		gamma_div[p, r, ] <- diversity_3ways(m = neighbs, flavor = 'gamma', 
											 dotd=T, dopd=T, dofd=F, abundance=T,
											 pddist = fiadist, fdmat = try_noproblem,
											 nnull = nnull,
											 phylo_spp = pnwphylo$tip.label, func_problem_spp = problemspp)
	}

}

close(pb)

cnames <- c('richness', 'shannon', 'evenness',
            'MPD_pa_z', 'MNTD_pa_z',
            'MPD_z', 'MNTD_z',
            'FRic', 'FEve', 'FDiv', 'FDis',
            'FRic_pa', 'FEve_pa', 'FDiv_pa', 'FDis_pa')

dimnames(gamma_div) <- list(cnames, paste(r,radii,sep='_'), NULL)

save(gamma_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/gamma_', slice, '.r'))