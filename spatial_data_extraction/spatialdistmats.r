# Spatial distances of all coordinate pairs.

bbs_coords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')

library(sp)
bbs_distmat <- spDists(x = with(bbs_coords, cbind(lon, lat)), longlat = TRUE)

# for interest, output number of pairwise comparisons that fall under each radius.
for (i in 100 * 1:10) {
	print(sum(bbs_distmat <= i))
}

# There are 10e6 entries in matrix, but only 2.7e6 within 1000 km.

# FIA

fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)

library(dplyr)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])
			
			
fia_distmat <- spDists(x = with(fiacoords, cbind(lon, lat)), longlat = TRUE)

# for interest, output number of pairwise comparisons that fall under each radius.
for (i in 100 * 1:10) {
	print(sum(fia_distmat <= i))
}

dimnames(bbs_distmat)[[1]] <- dimnames(bbs_distmat)[[2]] <- bbs_coords$rteNo
dimnames(fia_distmat)[[1]] <- dimnames(fia_distmat)[[2]] <- fiacoords$PLT_CN

save(bbs_distmat, file = '/mnt/research/nasabio/data/bbs/bbsspatialdistmat.r')
save(fia_distmat, file = '/mnt/research/nasabio/data/fia/fiaspatialdistmat.r')
