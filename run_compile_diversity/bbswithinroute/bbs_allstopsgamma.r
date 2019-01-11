# Calculate gamma diversity within route for BBS.

# Edited 09 Jan 2019: add 1 km
# Edited 11 Jan 2019: do in parallel

load('/mnt/research/nasabio/data/bbs/bbsworkspace_bystop_20072016.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(foreach)
library(doParallel)
library(sp)
library(vegan)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 999

radii <- c(1, 5, 10, 20)
stop_bounds <- rbind(c(25,26), c(20,31), c(13,37), c(1,50)) # First and last stops to be used for each radius (approx. 1,5,10,20 km)
route_ids <- unique(bbscov_oneyear$rteNo)

registerDoParallel(cores = length(radii))

# Loop through the radii and the unique route IDs.
gamma_div <- foreach (r = 1:length(radii)) %dopar% {
	gamma_div_r <- matrix(NA, nrow = length(route_ids), ncol = 11)
	for (p in 1:length(route_ids)) {
		
		# subset the stops within the radius along the route.
		rowidx <- bbscov_oneyear$rteNo == route_ids[p] & bbscov_oneyear$Stop >= stop_bounds[r, 1] & bbscov_oneyear$Stop <= stop_bounds[r, 2]
	  
		neighbs <- bbsmat_oneyear[rowidx, , drop = FALSE]
		gamma_div_r[p, ] <- diversity_3ways(m = neighbs, flavor = 'gamma', 
											 dotd=T, dopd=T, dofd=T, abundance=F,
											 pddist = ericdist, fddist = birdtraitdist,
											 nnull = nnull,
											 phylo_spp = NULL, func_problem_spp = NULL)

	}

	cnames <- c('richness', 'shannon', 'evenness',
				'MPD_pa_z', 'MNTD_pa_z',
				'MPD_z', 'MNTD_z',
				'MPDfunc_pa_z', 'MNTDfunc_pa_z',
				'MPDfunc_z', 'MNTDfunc_z')

	dimnames(gamma_div_r)[[2]] <- cnames
	gamma_div_r
}

gamma_div <- do.call(rbind, gamma_div)
gamma_div <- data.frame(rteNo = route_ids, radius = rep(radii, each = length(route_ids)), gamma_div)
gamma_div <- gamma_div[, apply(gamma_div, 2, function(x) any(!is.na(x)))]
write.csv(gamma_div, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/withinroute/bbs_withinroute_gamma.csv', row.names = FALSE)
