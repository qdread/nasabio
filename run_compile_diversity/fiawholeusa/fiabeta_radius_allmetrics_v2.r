# Calculate median beta-diversity of fia by radius.
# Modified 25 June to include PD and FD.
# Modified 09 Sep to load a different object
# Modified 04 Dec for the unfuzzed coordinates
# Modified 05 Dec to use the mean rather than the median
# Modified 06 Dec to use arcsin sqrt
# Modified 08 Jan 2018: for taxonomic only, and load only the slice that's needed at the time (parallel job).
# New version created 01 Mar 2018: keep only natural plots, don't use plantation plots
# Another version created 04 Apr 2018: do all metrics separately now that TD, PD, and FD are all done. 
# Yet another version created 28 Nov 2018: new version with new FIA data (macroplots removed) and slurm IDs
# Modified 14 Dec 2018: now uses median, not transformed mean.
# Modified 18 Dec 2018: uses the huge single sparse matrix.
# Yet another version created 02 Jan 2019: loads each plot and runs code on it

# Load FIA coordinates
load('/mnt/home/qdr/data/fiaworkspace_spatial_wholeusa_2018.r')

###############################################
# subset to keep only the natural plots
library(dplyr)
plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
fiacoords <- left_join(fiacoords, plantation)

# Edit 18 Dec: get rid of the non plantation rows now.
fiacoords <- filter(fiacoords, !plantation)

###############################################

# For each year and route number, get the median beta diversity within each radius.
# The pairwise table was constructed to only go up to 100 km, so that is all we will have.
radii <- c(5, 10, 20, 50, 75, 100) # in km
n_plot <- 119177

library(sp)
library(foreach)
library(doParallel)
registerDoParallel(cores = 24)

fp <- '/mnt/research/nasabio/data/fia/diversity/usa2018'

allmetrics2sparse <- function(x) {
	notna <- apply(x, 1, function(z) any(!is.na(z)))
	cbind(which(notna), x[notna, , drop = FALSE])
}


# Function 
neighbordivfromfullmatrix <- function(i) {
	lon <- fiacoords$lon[i]
	lat <- fiacoords$lat[i]
	load(file.path(fp, paste0('beta_', i, '.r')))
	beta_div <- allmetrics2sparse(beta_div)
	neighbors <- fiacoords[beta_div[,1], ]
	neighbordists <- spDistsN1(pts = cbind(neighbors$lon, neighbors$lat), pt = c(lon, lat), longlat = TRUE)
	commdat <- list()
	for (r in 1:length(radii)) {
		beta_incircle <- beta_div[neighbordists <= radii[r] ,-(1:2), drop = FALSE]
		beta_median <- apply(beta_incircle, 2, function(z) median(z[is.finite(z)]))

		commdat[[r]] <- c(PLT_CN = fiacoords$PLT_CN[i], radius = radii[r], beta_median)
						  
	}
	as.data.frame(do.call('rbind', commdat))
}

fia_beta <- foreach(i = 1:n_plot) %dopar% {
	if (i %% 1000 == 0) print(i)
	neighbordivfromfullmatrix(i)
}

fia_beta <- bind_rows(fia_beta)

write.csv(fia_beta, file = '/mnt/research/nasabio/data/fia/biodiversity_CSVs/updated_nov2018/fiausa_natural_beta.csv', row.names = FALSE)