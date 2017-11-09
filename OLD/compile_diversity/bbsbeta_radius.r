# Calculate median beta-diversity of bbs by radius.
# Modified 08 June to include PD and FD.

# Load bbs beta diversity and route coordinates
load('/mnt/research/nasabio/data/bbs/bbs_betadivtdpdfd_arraylist.r')
load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')

library(dplyr)

# For each year and route number, get the median beta diversity within each radius.
# The pairwise table was constructed to only go up to 300 km, so that is all we will have.
radii <- c(50, 75, 100, 150, 200, 300) # in km
years <- 1997:2016

# Find the right matrix in the lookup table, and for each plot, get the median pairwise beta-diversity within each radius.

library(sp)

neighbordivfromtable <- function(x) {
	neighbormat <- subset(bbscov, year == x$year)
	focalpointindex <- which(neighbormat$year == x$year & neighbormat$rteNo == x$rteNo)
	neighbordists <- spDistsN1(pts = cbind(neighbormat$lon, neighbormat$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- bbs_betadiv_arraylist[[which(years == x$year)]][neighbordists <= radii[i], focalpointindex, , drop = FALSE]
		commdat[[i]] <- c(radius = radii[i], apply(neighbors_incircle, 3, function(x) median(x[is.finite(x)])))

	}
	as.data.frame(do.call('rbind', commdat))
}


bbs_beta <- bbscov %>%
	filter(year >= 1997) %>%
	group_by(year, rteNo) %>%
	do(neighbordivfromtable(.))

bbs_beta <- bbs_beta %>% ungroup %>% left_join(bbscov) %>% select(c(1:2, 21:24, 3:20))
	
write.csv(bbs_beta, file = '/mnt/research/nasabio/data/bbs/bbs_beta_byroute.csv', row.names = FALSE)	
	