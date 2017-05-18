# Calculate gamma diversity within different radii for BBS.
# in contrast to FIA, we can only use groups within the same year. 1997-present.

# Split into 250 groups.
# Use precalculated matrix.
# Loop through each BBS plot.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 500 km of the target plot.
# Do median alpha, and total gamma, for all the neighbors within each radius.
# All others outside that radius get NA.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
source('~/code/fia/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.
load('/mnt/research/nasabio/data/bbs/birdtraitmat_clean.r')

library(sp)
library(vegan)
source('~/code/fia/fixpicante.r')

nnull <- 99

# run through all radii in each task, but split by year and also slice it up some.
radii <- c(50, 75, 100, 150, 200, 250, 300, 400, 500)
years <- 1997:2015 # 19 years, so we can have 13 slices per year. 19*13 = 247
n_slices <- 13
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

combos <- expand.grid(years, 1:n_slices)
year <- combos[task, 1]
slice <- combos[task, 2]

# Split the master bbs matrix by year.
bbs_year_mat <- fixedbbsmat_byroute[bbscov$year == year, ]
bbs_year_coords <- bbscov[bbscov$year == year, ]

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(bbs_year_mat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Declare structures to hold data

pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
gamma_div <- array(NA, dim = c(nrow(bbs_year_mat), length(radii), 17))

for (p in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, p)
	# Distances between target plot and all other plots.
  dist_p <- spDistsN1(pts=with(bbs_year_coords, cbind(lon, lat)), pt = c(bbs_year_coords$lon[p], bbs_year_coords$lat[p]), longlat = TRUE)
  
	for (r in 1:length(radii)) {
		neighbs <- bbs_year_mat[dist_p <= radii[r], ]
		gamma_div[p, r, ] <- diversity_3ways(m = neighbs, flavor = 'gamma', 
											 dotd=T, dopd=T, dofd=F, abundance=F,
											 pddist = ericdist, fdmat = birdtraitclean,
											 nnull = nnull,
											 phylo_spp = NULL, func_problem_spp = NULL)
	}

}

close(pb)

cnames <- c('richness', 'shannon', 'evenness',
            'MPD_pa_z', 'MNTD_pa_z',
            'MPD_z', 'MNTD_z',
            'FRic', 'FEve', 'FDiv', 'FDis',
            'FRic_pa', 'FEve_pa', 'FDiv_pa', 'FDis_pa')

dimnames(gamma_div) <- list(cnames, paste(r,radii,sep='_'), NULL)

save(gamma_div, file = paste0('/mnt/research/nasabio/data/bbs/diversity/gamma_', year, '_', slice, '.r'))