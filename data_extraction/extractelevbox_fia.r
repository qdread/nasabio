n_slices <- 5000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# FIA lat long coordinates
library(dplyr)
library(raster)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')

# Load precalculated distances
load('/mnt/research/nasabio/data/fia/distlogical.r')

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

stats_by_point <- list()

pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)
		   
for (i in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, i)
	extractBox(coords = with(fiacoords, cbind(lon, lat))[i,,drop=FALSE],
		   raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt',
		   radius = 300,
		   fp = '/mnt/research/nasabio/temp',
		   filetags = paste('fiatest', row.names(fiacoords)[i], sep = '_'))
	file_i <- paste0('/mnt/research/nasabio/temp/bbox_fiatest_', row.names(fiacoords)[i], '.tif')
	stats_by_point[[length(stats_by_point) + 1]] <- statsByRadius(file_i)
	if (file.exists(file_i)) deleteBox(file_i)
}	

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats/newslice/stats_',slice,'.r'))