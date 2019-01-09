# Calculate gamma diversity within route for BBS.

# Edited 09 Jan 2019: add 1 km

load('/mnt/research/nasabio/data/bbs/bbsworkspace_bystop_20072016.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(sp)
library(vegan)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 999

radii <- c(1, 5, 10, 20)
stop_bounds <- rbind(c(25,26), c(20,31), c(13,37), c(1,50)) # First and last stops to be used for each radius (approx. 1,5,10,20 km)
route_ids <- unique(bbscov_oneyear$rteNo)

# Declare structure to hold data
gamma_div <- replicate(length(radii), matrix(NA, nrow = length(route_ids), ncol = 11), simplify = FALSE)

# Loop through the three radii (5,10,20) and the unique route IDs.
for (r in 1:length(radii)) {
	cat('\nRadius: ', radii[r],'\n')
	pb <- txtProgressBar(1, length(route_ids), style = 3)
	for (p in 1:length(route_ids)) {
		setTxtProgressBar(pb, p)
		
		# subset the stops within the radius along the route.
		rowidx <- bbscov_oneyear$rteNo == route_ids[p] & bbscov_oneyear$Stop >= stop_bounds[r, 1] & bbscov_oneyear$Stop <= stop_bounds[r, 2]
	  
		neighbs <- bbsmat_oneyear[rowidx, , drop = FALSE]
		gamma_div[[r]][p, ] <- diversity_3ways(m = neighbs, flavor = 'gamma', 
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

	dimnames(gamma_div[[r]])[[2]] <- cnames
}

gamma_div <- do.call(rbind, gamma_div)
gamma_div <- data.frame(rteNo = route_ids, radius = radii, each = length(route_ids)), gamma_div)
write.csv(gamma_div, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/withinroute/bbs_withinroute_gamma.csv', row.names = FALSE)
