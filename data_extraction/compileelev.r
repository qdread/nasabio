# Compile elevation stats for the west coast into one file.
# modified 1 May: all, not just west coast.

elev_list <- list()

for (i in 1:100) {
	load(paste0('/mnt/research/nasabio/data/fia/elevstats/30m/stats_', i, '.r'))
	elev_list <- c(elev_list, stats_by_point)
}


library(dplyr)

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

fiacoords <- as.data.frame(fiacoords)			
			
for (i in 1:length(elev_list)) {
	elev_list[[i]] <- data.frame(fiacoords[i,], elev_list[[i]])
}

fia_elev_stats <- do.call('rbind', elev_list)
write.csv(fia_elev_stats, file = '/mnt/research/nasabio/data/dem/fia_elev_stats_noalaska.csv', row.names = FALSE)

# Added 3 May: also extract bbs

elev_list <- list()

for (i in 1:150) {
	load(paste0('/mnt/research/nasabio/data/bbs/elevstats/30m/stats_', i, '.r'))
	elev_list <- c(elev_list, stats_by_point)
}


library(dplyr)

fp <- '/mnt/research/nasabio'
bbscoords <- read.csv(file.path(fp, 'data/bbs/bbs_wgs84_coords_byroute.csv'), stringsAsFactors = FALSE)
bbsaea <- read.csv(file.path(fp, 'data/bbs/bbs_aea_coords_byroute.csv'), stringsAsFactors = FALSE)
load(file.path(fp, 'data/bbs/bbsworkspace_byroute.r'))

# paste all coordinates, year, and route number together.
bbsall <- data.frame(year=rep(NA,nrow(bbscoords)), rteNo=NA, x_aea=NA, y_aea=NA)
bbsall[which(!is.na(bbscoords[,1])), ] <- bbscov
bbsall <- transform(bbsall, rteNo=as.numeric(rteNo))
bbsall <- cbind(bbsall, bbscoords)

for (i in 1:length(elev_list)) {
	elev_list[[i]] <- data.frame(bbsall[i,], elev_list[[i]])
}

bbs_elev_stats <- do.call('rbind', elev_list)
write.csv(bbs_elev_stats, file = file.path(fp, 'data/dem/bbs_elev_stats.csv'), row.names = FALSE)