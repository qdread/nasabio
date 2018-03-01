# Calculate median beta-diversity of fia by radius.
# Modified 25 June to include PD and FD.
# Modified 09 Sep to load a different object
# Modified 04 Dec for the unfuzzed coordinates
# Modified 05 Dec to use the mean rather than the median
# Modified 06 Dec to use arcsin sqrt
# Modified 08 Jan 2018: for taxonomic only, and load only the slice that's needed at the time (parallel job).
# New version created 01 Mar 2018: keep only natural plots, don't use plantation plots

# Load FIA coordinates

slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

source('/mnt/research/nasabio/code/loadfiaall.r')

###############################################
# subset to keep only the natural plots
library(dplyr)
plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
fiacoords <- left_join(fiacoords, plantation)

###############################################

# Determine row indices for the slice of the matrix to be used.
n_slices <- 250
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
fiacoords$n <- 1:nrow(fiacoords)

# For each year and route number, get the median beta diversity within each radius.
# The pairwise table was constructed to only go up to 300 km, so that is all we will have.
radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300) # in km

library(sp)

neighbordivfromtable <- function(x) {
	neighbordists <- spDistsN1(pts = cbind(fiacoords$lon, fiacoords$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- beta_div_list[[x$n - rowidxmin + 1]][neighbordists <= radii[i] & !fiacoords$plantation, , drop = FALSE]
		commdat[[i]] <- c(radius = radii[i], 
						  apply(neighbors_incircle[, prop_vars, drop = FALSE], 2, function(x) sin(mean(asin(sqrt(x[is.finite(x)]))))^2),
						  apply(neighbors_incircle[, -(prop_vars), drop = FALSE], 2, function(x) mean(x[is.finite(x)])))
	}
	as.data.frame(do.call('rbind', commdat))
}


print(slice)
# Load the correct slice.
load(paste0('/mnt/research/nasabio/data/fia/diversity/usa/tdbeta_', slice, '.r'))

rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Identify which means need to be done on proportion transform.
prop_vars <- which(dimnames(beta_div_list[[1]])[[2]] %in% c('beta_td_pairwise', 'beta_td_sorensen', 'beta_td_pairwise_pa', 'beta_td_sorensen_pa'))


# Find the right matrix in the lookup table, and for each plot, get the median pairwise beta-diversity within each radius.

fia_beta <- fiacoords[rowidxmin:rowidxmax, ] %>%
	rowwise() %>%
	do(neighbordivfromtable(.))

fia_beta <- cbind(as.data.frame(fiacoords[rep(rowidxmin:rowidxmax, each = length(radii)), 'PLT_CN', drop = FALSE]), as.data.frame(fia_beta[,1:5]))

write.csv(fia_beta, file = paste0('/mnt/research/nasabio/data/fia/diversity/filtered/tdbeta_averages_', slice, '.csv'), row.names = FALSE)

#################################
# Compile the written csvs into one.
fia_beta_list <- list()

for (i in 1:250) {
	fia_beta_list[[i]] <- read.csv(paste0('/mnt/research/nasabio/data/fia/diversity/filtered/tdbeta_averages_', i, '.csv'), stringsAsFactors = FALSE)
}

fia_beta_list <- do.call('rbind', fia_beta_list)

plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
fia_beta_list <- subset(fia_beta_list, PLT_CN %in% plantation$PLT_CN[!plantation$plantation])

write.csv(fia_beta_list, file = '/mnt/research/nasabio/data/fia/fiausa_natural_betatd.csv', row.names = FALSE)