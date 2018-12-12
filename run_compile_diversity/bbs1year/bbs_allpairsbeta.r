# in each slice, load calculated bbs matrix.
# do a certain number of the pairwise comparisons. 

# Split into 2k groups
# Use precalculated matrix.
# Loop through each BBS route.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 500 km of the target plot.
# Do all pairwise taxonomic beta diversity between that plot and all its neighbors.
# All others outside that radius get NA.

# Modified 12 Dec 2018: edit maximum radius, number of iterations, and don't do pairwise things twice!

load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(sp)
library(vegan)
library(vegetarian, lib.loc = '/mnt/home/qdr/R/x86_64-pc-linux-gnu-library/3.5')
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 999 # edited 12dec2018


max_radius <- 100 # edited 12dec2018
n_slices <- 2000
slice <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + as.numeric(Sys.getenv('N1000')) * 1000

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(bbsmat_byroute_oneyear),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Initialize structures to hold data

# Only use progress bar if >1 plot
if (rowidxmin < rowidxmax) pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
beta_div <- list()

for (p in rowidxmin:rowidxmax) {
	if (rowidxmin < rowidxmax) setTxtProgressBar(pb, p)
	# Distances between target plot and all other plots.
	dist_p <- spDistsN1(pts=with(bbscov_oneyear, cbind(lon, lat)), pt = c(bbscov_oneyear$lon[p], bbscov_oneyear$lat[p]), longlat = TRUE)
	beta_div_p <- matrix(NA, nrow=nrow(bbsmat_byroute_oneyear), ncol=21)
		
	for (p2 in p:nrow(bbscov_oneyear)) {
		# Loop through all other FIA plots, check if plot is in radius
		# If plot is within radius, calculate diversity between that plot and target plot. 
		if (dist_p[p2] > 0 & dist_p[p2] <= max_radius) {
			beta_div_p[p2,] <- singlepair_beta(p1 = bbsmat_byroute_oneyear[p,], p2 = bbsmat_byroute_oneyear[p2,], 
												td=T, pd=T, fd=T, abundance=F,
												pddist=ericdist, fddist=birdtraitdist,
												nnull = nnull,
												phylo_spp = NULL, func_problem_spp = NULL)
		}
	}
	
	dimnames(beta_div_p)[[2]] <- c('beta_td_pairwise_pa', 'beta_td_sorensen_pa',
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
	beta_div[[length(beta_div) + 1]] <- beta_div_p
	
}

save(beta_div, file = paste0('/mnt/research/nasabio/data/bbs/diversity1year/beta_', slice, '.r'))
     
if (rowidxmin < rowidxmax) close(pb)

