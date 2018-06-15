# BBS beta-diversity within routes

# Edited 15 June 2018: Correct indexing error.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_bystop_20072016.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(sp)
library(vegan)
library(vegetarian)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 99 # Reduce to save time

n_slices <- 1000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

stop_bounds <- rbind(c(20,31), c(13,37), c(1,50)) # First and last stops to be used for each radius (approx. 5,10,20 km)
route_ids <- unique(bbscov_oneyear$rteNo)

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,length(route_ids),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

div_names <- c('beta_td_pairwise_pa', 'beta_td_sorensen_pa',
				'beta_td_pairwise', 'beta_td_sorensen',
				'beta_td_shannon', 
				'beta_pd_pairwise_pa', 'beta_pd_pairwise_pa_z',
				'beta_pd_nt_pa', 'beta_pd_nt_pa_z',
				'beta_pd_pairwise', 'beta_pd_pairwise_z',
				'beta_pd_nt', 'beta_pd_nt_z',
				'beta_fd_pairwise_pa', 'beta_fd_pairwise_pa_z',
				'beta_fd_nt_pa', 'beta_fd_nt_pa_z',
				'beta_fd_pairwise', 'beta_fd_pairwise_z',
				'beta_fd_nt', 'beta_fd_nt_z')
				
# Identify which means need to be done on proportion transform.
prop_vars <- which(div_names %in% c('beta_td_pairwise', 'beta_td_sorensen', 'beta_td_pairwise_pa', 'beta_td_sorensen_pa'))

radii <- c(5, 10, 20)

# Get pairwise metrics between: stops 1-49, and all stops with a greater number than that.
# 3 nested loops: rte = route, p = first plot, 1:49, p2 = second plot, (p+1):50

# Take the averages already within this script, so that we don't need to save the individual pairwise distances.

# Only use progress bar if >1 plot
if (rowidxmin < rowidxmax) pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
beta_div <- list()

for (rte in rowidxmin:rowidxmax) {
	if (rowidxmin < rowidxmax) setTxtProgressBar(pb, rte)
	# Distances between target plot and all other plots.
	beta_div_rte <- array(NA, dim=c(50,50,21))
	
	# Get indexes of rows that are in the focal route.
	rte_rows <- which(bbscov_oneyear$rteNo == route_ids[rte])
	
	for (p1 in 1:49) {
		for (p2 in (p1+1):50) {
			# Added 15 June: Catch errors and return NA
			beta_p1_p2 <- try(  singlepair_beta(p1 = bbsmat_oneyear[rte_rows[p1],], p2 = bbsmat_oneyear[rte_rows[p2],], 
												td=T, pd=T, fd=T, abundance=F,
												pddist=ericdist, fddist=birdtraitdist,
												nnull = nnull,
												phylo_spp = NULL, func_problem_spp = NULL), TRUE)
			beta_div_rte[p1, p2, ] <- if (!inherits(beta_p1_p2, 'try-error')) beta_p1_p2 else as.numeric(rep(NA, length(div_names)))
			
		}
	}
	
# Get mean beta-diversity of the point by applying a function to each 50x50 slice, along the 3rd dimension of beta_div_rte.	
	commdat <- list()
	for (r in 1:3) {
		neighbors_incircle <- beta_div_rte[stop_bounds[r,1]:stop_bounds[r,2], stop_bounds[r,1]:stop_bounds[r,2], , drop = FALSE]
		commdat[[r]] <- c(radius = radii[r], 
						  apply(neighbors_incircle[, , prop_vars, drop = FALSE], 3, function(x) sin(mean(asin(sqrt(x[is.finite(x)]))))^2),
						  apply(neighbors_incircle[, , -(prop_vars), drop = FALSE], 3, function(x) mean(x[is.finite(x)])))
	}
	commdat <- as.data.frame(do.call('rbind', commdat))
	names(commdat) <- c('radius', div_names)
	beta_div[[length(beta_div) + 1]] <- commdat
	
}

save(beta_div, file = paste0('/mnt/research/nasabio/data/bbs/diversitywithinroute/beta_', slice, '.r'))
     
if (rowidxmin < rowidxmax) close(pb)

