# Calculate median beta-diversity of fia by radius.
# Modified 25 June to include PD and FD.
# Modified 09 Sep to load a different object
# Modified 04 Dec for the unfuzzed coordinates
# Modified 05 Dec to use the mean rather than the median
# Modified 06 Dec to use arcsin sqrt

# Load bbs beta diversity and route coordinates
load('/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_listbymetric.r') # v. large (13gb on hard disk)
library(abind)
fia_betadiv_array <- abind(fia_betadiv_list, along = 3)
rm(fia_betadiv_list)

load('/mnt/research/nasabio/data/fia/fiaworkspace_nospatial.r')
source('/mnt/research/nasabio/code/loadfia.r')

library(dplyr)

# For each year and route number, get the median beta diversity within each radius.
# The pairwise table was constructed to only go up to 300 km, so that is all we will have.
radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300) # in km

# Identify which means need to be done on proportion transform.
prop_vars <- which(dimnames(fia_betadiv_array)[[3]] %in% c('beta_td_pairwise', 'beta_td_sorensen', 'beta_td_pairwise_pa', 'beta_td_sorensen_pa'))


# Find the right matrix in the lookup table, and for each plot, get the median pairwise beta-diversity within each radius.

library(sp)

neighbordivfromtable <- function(x) {
	focalpointindex <- which(fiacoords$PLT_CN == x$PLT_CN)
	neighbordists <- spDistsN1(pts = cbind(fiacoords$lon, fiacoords$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- fia_betadiv_array[neighbordists <= radii[i], focalpointindex, , drop = FALSE]
		commdat[[i]] <- c(radius = radii[i], 
						  apply(neighbors_incircle[, , prop_vars, drop = FALSE], 3, function(x) sin(mean(asin(sqrt(x[is.finite(x)]))))^2),
						  apply(neighbors_incircle[, , -(prop_vars), drop = FALSE], 3, function(x) mean(x[is.finite(x)])))
	}
	as.data.frame(do.call('rbind', commdat))
}


fia_beta <- fiacoords %>%
	rowwise() %>%
	do(neighbordivfromtable(.))

fia_beta <- cbind(as.data.frame(fiacoords[rep(1:nrow(fiacoords), each = length(radii)), c('PLT_CN', 'STATECD', 'COUNTYCD', 'PLOT')]), as.data.frame(fia_beta))
write.csv(fia_beta, file = '/mnt/research/nasabio/data/fia/fia_beta.csv', row.names = FALSE)