# Compile BBS diversity as a single value for each route, all years 2007-2016 pooled together
# QDR 26 Sep 2017 nasabioxgeo

########################################################################
# Alpha diversity

bbs_alphadiv <- read.csv('/mnt/research/nasabio/data/bbs/biodiversity_CSVs/bbs_alphadiv.csv', stringsAsFactors = FALSE)

library(dplyr)

# Replace coordinates with correct ones
bbscoords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv')

bbs_alphadiv <- bbscoords %>% left_join(bbs_alphadiv[,!names(bbs_alphadiv) %in% c('lon','lat','lon_aea','lat_aea')])
bbs_alphadiv <- filter(bbs_alphadiv, !is.na(richness))

# For each year and route number, get the median alpha diversity within each radius.
radii <- c(50, 75, 100, 150, 200, 300, 400, 500) # in km
radii <- c(50, 75, 100) # Added 13 Dec 2018 because newer version does not need to go above 100 km

library(sp)

neighbordiv <- function(x, dat) {
	neighbordists <- spDistsN1(pts = cbind(dat$lon, dat$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- dat[neighbordists <= radii[i], ]
		commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], richness = median(richness), MPD_pa_z = median(MPD_pa_z, na.rm=TRUE), MNTD_pa_z = median(MNTD_pa_z, na.rm=TRUE), MPDfunc_pa_z = median(MPDfunc_pa_z, na.rm = TRUE), MNTDfunc_pa_z = median(MNTDfunc_pa_z[is.finite(MNTDfunc_pa_z)])))
	}
	as.data.frame(do.call('rbind', commdat))
}

bbs_alpha <- bbs_alphadiv %>% 
	rowwise %>%
	do(neighbordiv(., dat = bbs_alphadiv))
	
bbs_alpha <- cbind(bbs_alphadiv[rep(1:nrow(bbs_alphadiv), each=length(radii)),c('rteNo','lon','lat','lon_aea','lat_aea')], bbs_alpha)
bbs_alpha <- arrange(bbs_alpha, rteNo, radius)	
write.csv(bbs_alpha, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/bbs_alpha_1year.csv', row.names = FALSE)	

########################################################################
# Gamma diversity
radii <- c(50, 75, 100, 150, 200, 300, 400, 500) # in km
radii <- c(50, 75, 100) # Added 13 Dec 2018 because newer version does not need to go above 100 km

bbs_gammadiv <- list()

for (i in radii) {
	load(paste0('/mnt/research/nasabio/data/bbs/diversity1year/gamma',i,'.r'))
	gamma_div <- cbind(radius = i, gamma_div)
	bbs_gammadiv[[length(bbs_gammadiv) + 1]] <- gamma_div
}

load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')

bbs_gamma <- cbind(bbscov_oneyear[rep(1:nrow(bbscov_oneyear), times=length(bbs_gammadiv)),c('rteNo','lon','lat','lon_aea','lat_aea')], do.call('rbind', bbs_gammadiv))

write.csv(bbs_gamma, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/bbs_gamma_1year.csv', row.names = FALSE)

########################################################################
# Beta diversity older method

# Make an array of route x route x metric.

bbs_betadiv <- list()

for (i in 1:2000) {
	load(paste0('/mnt/research/nasabio/data/bbs/diversity1year/beta_', i, '.r'))
	bbs_betadiv[[i]] <- beta_div
	print(i)
}

bbs_betadiv <- do.call('c', bbs_betadiv) # Flatten into one list. # All are 3089 rows x 21 columns.

# Convert this to a 3089 x 3089 x 21 array.
library(abind)
bbs_betadiv_array <- abind(bbs_betadiv, along = 0)

# Copy the mirror image of the upper triangle into the lower
copy_upper <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
}

for (i in 1:dim(bbs_betadiv_array)[3]) bbs_betadiv_array[, , i] <- copy_upper(bbs_betadiv_array[, , i])

save(bbs_betadiv_array, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/bbs_betadivtdpdfd_array_1year.r')

########################################################################
# Beta diversity newer method

# Compile "final" beta-partition metrics for bbs. (Baselga only, Sorensen only, taxonomic and phylogenetic)


bbs_betapart_route <- list()

radii <- c(50, 75, 100, 150, 200, 300, 400, 500) * 1000
slices_per_radius <- 250

for (r in 1:length(radii)) {
	list_r <- list()
	
	for (slice in 1:slices_per_radius) {
		dat <- read.csv(paste0('/mnt/research/nasabio/data/bbs/diversity1year/betapartfinal_',as.character(as.integer(radii[r])),'_',slice,'.csv'))
		list_r[[slice]] <- dat
	}

	bbs_betapart_route[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))
}

bbs_betapart_route <- do.call('rbind', bbs_betapart_route)
bbs_betapart_route <- cbind(bbs_alphadiv[rep(1:nrow(bbs_alphadiv), each=6),c('rteNo','lon','lat','lon_aea','lat_aea')], bbs_betapart_route)
write.csv(bbs_betapart_route, file = '/mnt/research/nasabio/data/bbs/bbs_betapart_1year.csv', row.names = FALSE)