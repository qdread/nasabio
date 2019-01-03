# in each slice, load calculated bbs matrix.
# do a certain number of the pairwise comparisons. 
# for BBS, this only has to be done within a year.

# Split into 250 groups.
# Use precalculated matrix.
# Loop through each BBS route.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 500 km of the target plot.
# Do all pairwise taxonomic beta diversity between that plot and all its neighbors.
# All others outside that radius get NA.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
source('~/code/fia/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(sp)
library(vegan)
library(vegetarian)
source('~/code/fia/fixpicante.r')

nnull <- 99 # Reduce to save time


max_radius <- 300e3 # 500 km
years <- 1997:2015
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

# Initialize structures to hold data

pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
beta_div <- list()

for (p in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, p)
	# Distances between target plot and all other plots.
	dist_p <- spDistsN1(pts=with(bbs_year_coords, cbind(lon, lat)), pt = c(bbs_year_coords$lon[p], bbs_year_coords$lat[p]), longlat = TRUE)
	beta_div_p <- matrix(NA, nrow=nrow(bbs_year_mat), ncol=19)
		
	for (p2 in 1:nrow(bbs_year_coords)) {
		# Loop through all other BBS routes, check if plot is in radius
		# If plot is within radius, calculate diversity between that plot and target plot. 
		if (dist_p[p2] > 0 & dist_p[p2] <= max_radius) {
			beta_div_p[p2,] <- singlepair_beta(p1 = bbs_year_mat[p,], p2 = bbs_year_mat[p2,], 
												td=T, pd=F, fd=F, abundance=F,
												pddist=ericdist, fddist=birdtraitdist,
												nnull = nnull,
												phylo_spp = NULL, func_problem_spp = NULL)
		}
	}
	
	dimnames(beta_div_p)[[2]] <- c('beta_td_pairwise_pa', 'beta_td_pairwise', 'beta_td_shannon', 
	                                'beta_pd_pairwise_pa', 'beta_pd_pairwise_pa_z',
	                                'beta_pd_nt_pa', 'beta_pd_nt_pa_z',
	                                'beta_pd_pairwise', 'beta_pd_pairwise_z',
	                                'beta_pd_nt', 'beta_pd_nt_pa_z',
	                                'beta_fd_pairwise_pa', 'beta_fd_pairwise_pa_z',
	                                'beta_fd_nt_pa', 'beta_fd_nt_pa_z',
	                                'beta_fd_pairwise', 'beta_fd_pairwise_z',
	                                'beta_fd_nt', 'beta_fd_nt_pa_z')
	beta_div[[length(beta_div) + 1]] <- beta_div_p
	
}

save(beta_div, file = paste0('/mnt/research/nasabio/data/bbs/diversity/tdbeta_', year,'_', slice, '.r'))
     
close(pb)

