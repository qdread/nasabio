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

bbscov <- bbsgrps %>%
	left_join(bbscoords %>% mutate(rteNo=as.character(rteNo))) %>%
	setNames(nm = c('year','rteNo','Stop','lon','lat','lon_aea','lat_aea'))

bbscov$Stop <- as.integer(gsub('Stop','',bbscov$Stop))
bbscov$rteNo <- as.integer(bbscov$rteNo)
# For optimization purposes, convert covariates to a matrix.
bbscovmat <- as.matrix(bbscov)

######
# Combine 2007-2016 into a single year.
# Combine each stop across years.

# Also do 2001-2011 in a single year, in case we want to distinguish between the two.

consolidate_years <- function(x) {
	mat_x <- fixedbbsmat[x$rowidx, , drop = FALSE]
	as.numeric(apply(mat_x, 2, sum) > 0)
}

bbs_consol <- bbscov %>%
	mutate(rowidx = 1:nrow(bbscov)) %>%
	filter(year >= 2007 & year <= 2016) %>%
	group_by(rteNo, Stop, lon, lat, lon_aea, lat_aea) %>%
	do(x = consolidate_years(.))
	
bbsmat_oneyear <- do.call('rbind', bbs_consol$x)
dimnames(bbsmat_oneyear)[[2]] <- dimnames(fixedbbsmat)[[2]]

bbscov_oneyear <- bbs_consol %>% select(-x)
bbscovmat_oneyear <- as.matrix(bbscov_oneyear)

save(bbscov_oneyear, bbscovmat_oneyear, bbsmat_oneyear, file = '/mnt/research/nasabio/data/bbs/bbsworkspace_bystop_20072016.r')


bbs_consol2001 <- bbscov %>%
	mutate(rowidx = 1:nrow(bbscov)) %>%
	filter(year >= 2001 & year <= 2011) %>%
	group_by(rteNo, Stop, lon, lat, lon_aea, lat_aea) %>%
	do(x = consolidate_years(.))
	
bbsmat_oneyear <- do.call('rbind', bbs_consol2001$x)
dimnames(bbsmat_oneyear)[[2]] <- dimnames(fixedbbsmat)[[2]]

bbscov_oneyear <- bbs_consol2001 %>% select(-x)
bbscovmat_oneyear <- as.matrix(bbscov_oneyear)

save(bbscov_oneyear, bbscovmat_oneyear, bbsmat_oneyear, file = '/mnt/research/nasabio/data/bbs/bbsworkspace_bystop_20012011.r')