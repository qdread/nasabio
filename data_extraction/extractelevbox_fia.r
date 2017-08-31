n_slices <- 5000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
raster_file_name <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem.vrt'
scratch_path <- paste0('/mnt/ffs17/users/', Sys.getenv('USER'))

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
		   
for (j in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, j)
	extractBox(coords = with(fiacoords, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = raster_file_name,
		   radius = 300,
		   fp = scratch_path,
		   filetags = paste('fiatest', row.names(fiacoords)[j], sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_fiatest_', row.names(fiacoords)[j], '.tif')
	stats_by_point[[length(stats_by_point) + 1]] <- statsByRadius(file_j)
	if (file.exists(file_j)) deleteBox(file_j)
}	

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats/newslice/stats_',slice,'.r'))
