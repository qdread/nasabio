# Calculate gamma diversity within different radii for BBS.
# in contrast to FIA, we can only use groups within the same year. 1997-present.

# Split into 250 groups.
# Use precalculated matrix.
# Loop through each BBS plot.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 100 km of the target plot.
# Do median alpha, and total gamma, for all the neighbors within each radius.
# All others outside that radius get NA.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

# Replace AOU codes in the trait matrix with actual species names.

library(sp)
library(vegan)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 999 # Edited 12 Dec 2018

# run through all radii in each task, but split by year and also slice it up some.
radii <- c(50, 75, 100) # edited 12 dec 2018
task <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
r <- radii[task]

# Declare structures to hold data

pb <- txtProgressBar(1, nrow(bbsmat_byroute_oneyear), style = 3)
gamma_div <- matrix(NA, nrow = nrow(bbsmat_byroute_oneyear), ncol = 11)

for (p in 1:nrow(bbsmat_byroute_oneyear) ) {
	setTxtProgressBar(pb, p)
	# Distances between target plot and all other plots.
  dist_p <- spDistsN1(pts=with(bbscov_oneyear, cbind(lon, lat)), pt = c(bbscov_oneyear$lon[p], bbscov_oneyear$lat[p]), longlat = TRUE)
  
		neighbs <- bbsmat_byroute_oneyear[dist_p <= r, , drop = FALSE]
		gamma_div[p, ] <- diversity_3ways(m = neighbs, flavor = 'gamma', 
											 dotd=T, dopd=T, dofd=T, abundance=F,
											 pddist = ericdist, fddist = birdtraitdist,
											 nnull = nnull,
											 phylo_spp = NULL, func_problem_spp = NULL)

}

close(pb)

cnames <- c('richness', 'shannon', 'evenness',
            'MPD_pa_z', 'MNTD_pa_z',
            'MPD_z', 'MNTD_z',
            'MPDfunc_pa_z', 'MNTDfunc_pa_z',
            'MPDfunc_z', 'MNTDfunc_z')

dimnames(gamma_div)[[2]] <- cnames

save(gamma_div, file = paste0('/mnt/research/nasabio/data/bbs/diversity1year/gamma',as.character(as.integer(r)), '.r'))
