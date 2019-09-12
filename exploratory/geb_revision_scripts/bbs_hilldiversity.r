# Code to only run beta-diversity in effective species numbers units for BBS
# Hopefully can be run very quickly.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')
source('/mnt/research/nasabio/code/pairwise_beta_focal.r')
load('/mnt/research/nasabio/data/bbs/bbspdfddist.r') # Phy and Func distance matrices.

library(sp)
library(vegetarian, lib.loc = '/mnt/home/qdr/R/x86_64-pc-linux-gnu-library/3.6')

max_radius <- 50

rowidxmin <- 1
rowidxmax <- nrow(bbsmat_byroute_oneyear)

# Initialize structures to hold data

# Only use progress bar if >1 plot
if (rowidxmin < rowidxmax) pb <- txtProgressBar(rowidxmin, rowidxmax, style = 3)
beta_div <- matrix(NA, nrow = rowidxmax, ncol = rowidxmax)

for (p in rowidxmin:rowidxmax) {
	if (rowidxmin < rowidxmax) setTxtProgressBar(pb, p)
	# Distances between target plot and all other plots.
	dist_p <- spDistsN1(pts=with(bbscov_oneyear, cbind(lon, lat)), pt = c(bbscov_oneyear$lon[p], bbscov_oneyear$lat[p]), longlat = TRUE)
		
	for (p2 in p:nrow(bbscov_oneyear)) {
		# If plot is within radius, calculate diversity between that plot and target plot. 
		if (dist_p[p2] > 0 & dist_p[p2] <= max_radius) {
			pair <- bbsmat_byroute_oneyear[c(p, p2), ]
			beta_div[p, p2] <- d(pair, lev = 'beta', q = 0)
		}
	}
}

save(beta_div, file = '/mnt/research/nasabio/data/bbs/diversity1year/additivebeta.RData')
     
if (rowidxmin < rowidxmax) close(pb)

