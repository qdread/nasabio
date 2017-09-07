# Compile summary statistics: elevation and bioclim.
# BBS and FIA
# 04 Sep 2017 QDR

##############################
# BBS

get_stats <- function(fp, prefix, n, suffix) {
	file_names <- paste0(prefix, n, suffix)
	all_stats <- list()
	for (i in file_names) {
		load(file.path(fp, i))
		all_stats[[length(all_stats) + 1]] <- stats_by_point
	}
	do.call('c', all_stats)
}

bbs_elev_stats <- get_stats('/mnt/research/nasabio/data/bbs/elevstats/30m', 'stats_', 1:1000, '.r')
bbs_bioclim5k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'bioclim5k_', 1:500, '.r')
bbs_bioclim1k_stats <- get_stats('/mnt/research/nasabio/data/bbs/climstats', 'bioclim1k_', 1:1000, '.r')

# Combine these, dealing with the fact that some list entries are NA.

bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv', stringsAsFactors = FALSE)
bbs_5k_missing <- bbs_bioclim5k_stats[[1]]
bbs_5k_missing[, 3:7] <- NA
bbs_1k_missing <- bbs_bioclim1k_stats[[1]]
bbs_1k_missing[, 3:7] <- NA

bbs_all_stats <- list()

for (i in 1:length(bbs_elev_stats)) {
	bbs_5k_i <- if (class(bbs_bioclim5k_stats[[i]]) == 'data.frame') bbs_bioclim5k_stats[[i]] else bbs_5k_missing
	bbs_5k_i$variable <- paste0('bio', 1:19)
	bbs_1k_i <- if (class(bbs_bioclim1k_stats[[i]]) == 'data.frame') bbs_bioclim1k_stats[[i]] else bbs_1k_missing
	bbs_1k_i$variable <- paste0('bio', 1:19,'_1k')
	bbs_elev_i <- data.frame(radius = bbs_elev_stats[[i]][,'radius'], variable = 'elevation', bbs_elev_stats[[i]][,c('mean','sd','min','max')], n = NA)
	bbs_all_stats[[i]] <- cbind(bbsll[i,], rbind(bbs_elev_i, bbs_5k_i, bbs_1k_i))
}

bbs_all_stats <- do.call('rbind', bbs_all_stats)

write.csv(bbs_all_stats, file = '/mnt/research/nasabio/data/bbs/bbs_geodiversity_stats.csv', row.names = FALSE)

##############################
# FIA

fia_elev_stats <- get_stats('/mnt/research/nasabio/data/fia/elevstats/newslice', 'stats_', 1:5000, '.r')
fia_bioclim5k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'bioclim5k_', 1:1000, '.r')
fia_bioclim1k_stats <- get_stats('/mnt/research/nasabio/data/fia/climstats', 'bioclim1k_', 1:5000, '.r')

# FIA lat long coordinates
library(dplyr)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

fia_5k_missing <- fia_bioclim5k_stats[[10001]]
fia_5k_missing[, 3:7] <- NA
fia_1k_missing <- fia_bioclim1k_stats[[10001]]
fia_1k_missing[, 3:7] <- NA
fia_elev_missing <- fia_elev_stats[[10001]]
fia_elev_missing[, 3:7] <- NA

fia_all_stats <- list()

for (i in 1:length(fia_elev_stats)) {
	fia_5k_i <- if (class(fia_bioclim5k_stats[[i]]) == 'data.frame') fia_bioclim5k_stats[[i]] else fia_5k_missing
	fia_5k_i$variable <- paste0('bio', 1:19)
	fia_1k_i <- if (class(fia_bioclim1k_stats[[i]]) == 'data.frame') fia_bioclim1k_stats[[i]] else fia_1k_missing
	fia_1k_i$variable <- paste0('bio', 1:19,'_1k')
	fia_elev_i <- if(class(fia_elev_stats[[i]]) == 'data.frame') fia_elev_stats[[i]] else fia_elev_missing
	fia_elev_i$variable <- 'elevation'
	fia_all_stats[[i]] <- cbind(as.data.frame(fiacoords[i,]), rbind(fia_elev_i, fia_5k_i, fia_1k_i))
}			

fia_all_stats <- do.call('rbind', fia_all_stats)

write.csv(fia_all_stats, file = '/mnt/research/nasabio/data/fia/fia_geodiversity_stats.csv', row.names = FALSE)
