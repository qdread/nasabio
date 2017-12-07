# in each slice, load calculated fia matrix.
# do a certain number of the pairwise comparisons. 
# at first, just do taxonomic.

# Split into 10k groups.
# Use precalculated matrix.
# Loop through each FIA plot.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 300 km of the target plot.
# Do all pairwise taxonomic beta diversity between that plot and all its neighbors.
# All others outside that radius get NA.

load('/mnt/research/nasabio/data/fia/fiaworkspace_nospatial.r')
source('/mnt/research/nasabio/code/loadfia.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')

library(sp)
library(vegan)
library(vegetarian)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 99 # Reduce to save time

trydist <- as.matrix(trydist)

max_radius <- 300 # Do 300 km for now.
n_slices <- 10000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaplotmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Initialize structures to hold data

pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
beta_div <- list()

for (p in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, p)
	# Distances between target plot and all other plots.
	dist_p <- spDistsN1(pts=with(fiacoords, cbind(lon, lat)), pt = c(fiacoords$lon[p], fiacoords$lat[p]), longlat = TRUE)
	beta_div_p <- matrix(NA, nrow=nrow(fiaplotmat), ncol=21)
		
	for (p2 in 1:nrow(fiaplotmat)) {
		# Loop through all other FIA plots, check if plot is in radius
		# If plot is within radius, calculate diversity between that plot and target plot. 
		if (dist_p[p2] > 0 & dist_p[p2] <= max_radius) {
			beta_div_p[p2,] <- singlepair_beta(p1 = fiaplotmat[p,], p2 = fiaplotmat[p2,], 
												td=T, pd=T, fd=T, abundance=T,
												pddist=fiadist, fddist=trydist,
												nnull = nnull,
												phylo_spp = pnwphylo$tip.label, func_problem_spp = NULL)
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

save(beta_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/unfuzzed/beta_', slice, '.r'))
     
close(pb)

