# Calculate median beta-diversity of bbs by radius. (old method)
# Modified 08 June to include PD and FD.
# Modified 05 Dec to use the mean rather than the median

# Load bbs beta diversity and route coordinates
load('/mnt/research/nasabio/data/bbs/bbs_betadivtdpdfd_array_1year.r')
load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')

library(dplyr)

# For each year and route number, get the median beta diversity within each radius.
# The pairwise table was constructed to only go up to 300 km, so that is all we will have.
radii <- c(50, 75, 100, 150, 200, 300) # in km

# Find the right matrix in the lookup table, and for each plot, get the median pairwise beta-diversity within each radius.

library(sp)

neighbordivfromtable <- function(x) {
	focalpointindex <- which(bbscov_oneyear$rteNo == x$rteNo)
	neighbordists <- spDistsN1(pts = cbind(bbscov_oneyear$lon, bbscov_oneyear$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- bbs_betadiv_array[neighbordists <= radii[i], focalpointindex, , drop = FALSE]
		commdat[[i]] <- c(radius = radii[i], apply(neighbors_incircle, 3, function(x) mean(x[is.finite(x)])))

	}
	as.data.frame(do.call('rbind', commdat))
}


bbs_beta <- bbscov_oneyear %>%
	rowwise %>%
	do(neighbordivfromtable(.))

bbs_beta <- cbind(bbscov_oneyear[rep(1:nrow(bbscov_oneyear), each=length(radii)),], bbs_beta)

write.csv(bbs_beta, file = '/mnt/research/nasabio/data/bbs/bbs_betatdpdfd_1year.csv', row.names = FALSE)	
	