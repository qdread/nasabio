n_slices <- 250
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
raster_file_name <- '/mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl.vrt'
scratch_path <- Sys.getenv('TMPDIR')

# FIA lat long coordinates
library(dplyr)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')

load('/mnt/research/nasabio/data/precalcdist/distlogical_fia_nightlight.r')

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(bbsll),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Radii of circles where we want summary stats
radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500) # in km

# Call extraction function, specifying the raster from which to extract data.
stats_by_point <- list()

pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)
		   
for (j in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, j)
	extractBox(coords = with(fiacoords, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = raster_file_name,
		   radius = 500,
		   fp = scratch_path,
		   filetags = paste('fianight', row.names(bbsll)[j], sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_fianight_', row.names(bbsll)[j], '.tif')
	stats_by_point[[length(stats_by_point) + 1]] <- statsByRadius(file_j, radii = radii, is_brick = FALSE)
	if (file.exists(file_j)) deleteBox(file_j)
}	

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/geostats/nightlight_',slice,'.r'))
