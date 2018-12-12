# FIA all pairs beta-diversity
# Edited 21 Dec for entire USA
# Edited 23 Dec: only do a single point per slice
# Edited 29 Dec (hpcc version only): max radius is now 200 to save time
# Edited 11 Jan 2018: Move all files to SCRATCH and/or TMPDIR.
# Edited 18 Jan 2018: Deal with jobs being over 100000.
# Edited 09 Feb 2018: Update scratch path
# Edited 21 Mar 2018: another update to scratch path
# Edited 26 Nov 2018: file paths update
# Edited 27 Nov 2018: get rid of parsing and update ID for SLURM
# Edited 12 Dec 2018: make four modifications: 1 - decrease_max radius, 2 - increase nnull, 3 - remove the plantation plots to save more time, 4 - don't do all pairwise twice!

# Use precalculated matrix.
# One FIA plot per slice.
# First, calculate spatial distances between that plot and all other plots.
# Identify the plots within 200 km of the target plot.
# Do all pairwise taxonomic beta diversity between that plot and all its neighbors.
# All others outside that radius get NA.

load('/mnt/ffs17/groups/nasabio/fiaworkspace_nospatial_wholeusa_2018.r')
load('/mnt/home/qdr/data/fiaworkspace_spatial_wholeusa_2018.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
source('/mnt/research/nasabio/code/nofuncspp.r')

# Added 12 Dec 2018: get rid of plantations
fiacoords <- fiacoords[!fiaplantation$plantation, ]
fiaplotmat <- fiaplotmat[!fiaplantation$plantation, ]

library(sp)
library(vegan)
library(vegetarian, lib.loc = '/mnt/home/qdr/R/x86_64-pc-linux-gnu-library/3.5')
source('/mnt/research/nasabio/code/fixpicante.r')

nnull <- 999 # edited 12dec2018

trydist <- as.matrix(trydist)

max_radius <- 100 # edited 12dec2018
p <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + as.numeric(Sys.getenv('N1000')) * 1000
print(p)

# Distances between target plot and all other plots.
dist_p <- spDistsN1(pts=with(fiacoords, cbind(lon, lat)), pt = c(fiacoords$lon[p], fiacoords$lat[p]), longlat = TRUE)
beta_div <- matrix(NA, nrow=nrow(fiaplotmat), ncol=21)

pb <- txtProgressBar(p, nrow(fiaplotmat), style = 3)
	
for (p2 in p:nrow(fiaplotmat)) {
	setTxtProgressBar(pb, p2)
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

close(pb)

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

save(beta_div, file = paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/beta_', as.integer(p), '.r'))
