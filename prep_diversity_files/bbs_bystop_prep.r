# Create input data for BBS at the stop level
# QDR, 17 Apr 2017
# Project: NASABioXGeo

# 1. Load community matrix.

bbsspp <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/specieslist.csv', stringsAsFactors = FALSE)
load('/mnt/research/aquaxterra/DATA/raw_data/BBS/bbsmat_bystop2016.r') # Load fixed bbsmat. This loads both byroute and bystop.
# Quick correction to fix two birds that aren't in the phylogeny. Just get rid of the eastern yellow wagtail since it's probably only in Alaska anyway.
fixedbbsmat[, which(sppids == 5739)] <- fixedbbsmat[, which(sppids == 5738)] + fixedbbsmat[, which(sppids == 5739)]
fixedbbsmat[, which(sppids == 5738)] <- 0
fixedbbsmat[, which(sppids == 6960)] <- 0

bbscoords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv')

birdtrait <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/birdtraitmerged.csv', stringsAsFactors = FALSE)
birdtrait[birdtrait == -999] <- NA

library(dplyr)

# Set dimnames of community matrix

dimnames_matrix <- rep(NA, length(sppids))
for (i in 1:length(sppids)) {
	names_i <- bbsspp[bbsspp$AOU == sppids[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_matrix[i] <- names_i[1]
}

dimnames(fixedbbsmat)[[2]] <- dimnames_matrix

# Remove nocturnal species, columns with zero sum, and routes that don't have a midpoint.
# Rows with zero sum will have to be removed later.
has_coords <- bbsgrps$rteNo %in% bbscoords$rteNo
ns <- colSums(fixedbbsmat)
nocturnalbirds <- birdtrait$Latin_Name[birdtrait$Nocturnal == 1]
fixedbbsmat <- fixedbbsmat[has_coords, !(dimnames(fixedbbsmat)[[2]] %in% nocturnalbirds) & ns != 0]

bbsgrps <- bbsgrps[has_coords, ]

#################### updated to this point.

# Calculate distances and identities of all neighbors within the maximum radius
names(bbsgrps_byroute) <- c('year','rteNo','lon','lat','lon_aea','lat_aea')
bbscov <- bbsgrps_byroute

# For optimization purposes, convert covariates to a matrix.
bbscovmat <- as.matrix(bbscov)


bbsnhb_list <- as.data.frame(bbscovmat) %>% group_by(year) %>% do(l = getNeighbors(., radius = 5e5)) 
# Flatten this into one list
bbsnhb_r <- do.call('c', bbsnhb_list$l)

save(bbsnhb_r, bbscov, bbscovmat, fixedbbsmat_byroute, file = '/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
write.csv(fixedbbsmat_byroute, file = '/mnt/research/nasabio/data/bbs/bbs_plot_matrix.csv', row.names = FALSE)

######
# Combine 2007-2016 into a single year.

consolidate_years <- function(x) {
	mat_x <- fixedbbsmat_byroute[x$rowidx, , drop = FALSE]
	as.numeric(apply(mat_x, 2, sum) > 0)
}

bbs_consol <- bbscov %>%
	mutate(rowidx = 1:nrow(bbscov)) %>%
	filter(year >= 2007) %>%
	group_by(rteNo, lon, lat, lon_aea, lat_aea) %>%
	do(x = consolidate_years(.))
	
bbsmat_byroute_oneyear <- do.call('rbind', bbs_consol$x)
dimnames(bbsmat_byroute_oneyear)[[2]] <- dimnames(fixedbbsmat_byroute)[[2]]

bbscov_oneyear <- bbs_consol %>% select(-x)
bbscovmat_oneyear <- as.matrix(bbscov_oneyear)

bbsnhb_list_oneyear <- getNeighbors(dat = as.data.frame(bbscovmat_oneyear), radius = 5e5)

save(bbsnhb_list_oneyear, bbscov_oneyear, bbscovmat_oneyear, bbsmat_byroute_oneyear, file = '/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')