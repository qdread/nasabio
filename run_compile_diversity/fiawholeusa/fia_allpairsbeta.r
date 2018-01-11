# FIA all pairs beta-diversity
# Edited 21 Dec for entire USA
# Edited 23 Dec: only do a single point per slice
# Edited 29 Dec (hpcc version only): max radius is now 200 to save time
# Edited 11 Jan 2018: Move all files to SCRATCH and/or TMPDIR.

# Use precalculated matrix.
# One FIA plot per slice.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 300 km of the target plot.
# Do all pairwise taxonomic beta diversity between that plot and all its neighbors.
# All others outside that radius get NA.

load('/mnt/ls15/scratch/groups/plz-lab/NASA/fiaworkspace_nospatial_wholeusa.r')
source('/mnt/research/nasabio/code/loadfiaall.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
source('/mnt/research/nasabio/code/nofuncspp.r')

library(sp)
library(vegan)
library(vegetarian)
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 99 # Reduce to save time

trydist <- as.matrix(trydist)

max_radius <- 200 
p <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# Distances between target plot and all other plots.
dist_p <- spDistsN1(pts=with(fiacoords, cbind(lon, lat)), pt = c(fiacoords$lon[p], fiacoords$lat[p]), longlat = TRUE)
beta_div <- matrix(NA, nrow=nrow(fiaplotmat), ncol=21)
	
for (p2 in 1:nrow(fiaplotmat)) {
	# Loop through all other FIA plots, check if plot is in radius
	# If plot is within radius, calculate diversity between that plot and target plot. 
	if (dist_p[p2] > 0 & dist_p[p2] <= max_radius) {
		beta_div[p2,] <- singlepair_beta(p1 = fiaplotmat[p,], p2 = fiaplotmat[p2,], 
											td = TRUE, pd = TRUE, fd = TRUE, abundance = TRUE,
											pddist=fiadist, fddist=trydist,
											nnull = nnull,
											phylo_spp = fullphylo$tip.label, func_problem_spp = nofuncspp)
	}
}

dimnames(beta_div)[[2]] <- c('beta_td_pairwise_pa', 'beta_td_sorensen_pa',
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

save(beta_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/usa/beta_', p, '.r'))
