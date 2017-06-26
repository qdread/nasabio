# Calculate median beta-diversity of fia by radius.
# Modified 25 June to include PD and FD.

# Load bbs beta diversity and route coordinates
load('/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_array.r') # v. large (13gb on hard disk)
load('/mnt/research/nasabio/data/fia/fiaworkspace2.r')

library(dplyr)

# For each year and route number, get the median beta diversity within each radius.
# The pairwise table was constructed to only go up to 300 km, so that is all we will have.
radii <- c(10, 20, 50, 75, 100, 150, 200, 300) # in km

# Find the right matrix in the lookup table, and for each plot, get the median pairwise beta-diversity within each radius.

library(sp)

neighbordivfromtable <- function(x) {
	focalpointindex <- which(fiacoords$PLT_CN == x$PLT_CN)
	neighbordists <- spDistsN1(pts = cbind(fiacoords$lon, fiacoords$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- fia_betadiv_array[neighbordists <= radii[i], focalpointindex, , drop = FALSE]
		commdat[[i]] <- c(radius = radii[i], apply(neighbors_incircle, 3, function(x) median(x[is.finite(x)])))

	}
	as.data.frame(do.call('rbind', commdat))
}


fia_beta <- fiacoords %>%
	rowwise() %>%
	do(neighbordivfromtable(.))

#fia_beta <- cbind(fiacoords, fia_beta)
	
#write.csv(fia_beta, file = '/mnt/research/nasabio/data/fia/fia_beta.csv', row.names = FALSE)	

load('/mnt/research/nasabio/data/fia/fiaworkspace2.r')	
load('/mnt/research/nasabio/data/fia/fia_betaobj.r')	

fia_beta <- cbind(fiacoords[rep(1:8, times = nrow(fiacoords)),], fia_beta)
write.csv(fia_beta, file = '/mnt/research/nasabio/data/fia/fia_beta.csv', row.names = FALSE)