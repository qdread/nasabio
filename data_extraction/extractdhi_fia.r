n_slices <- 250
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
layers <- c('fpar', 'gpp', 'lai8', 'ndvi')
raster_file_name <- paste0('/mnt/research/nasabio/data/dhi/dhi_', layers, '.vrt')
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

mat_list <- list()
for (ii in layers) {
	load(paste0('/mnt/research/nasabio/data/precalcdist/distlogical_fia_dhi',ii,'.r'))
	mat_list[[length(mat_list) + 1]] <- idx_r
}

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Radii of circles where we want summary stats
radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500) # in km

# Call extraction function, specifying the raster from which to extract data.
stats_by_point <- list()

pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)
		   
for (j in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, j)
	point_stat <- list()
	for (ii in 1:length(layers)) {
		idx_r <- mat_list[[ii]]
		extractBox(coords = with(fiacoords, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = raster_file_name[[ii]],
		   radius = 500,
		   fp = scratch_path,
		   filetags = paste('fiadhi', layers[[ii]], row.names(bbsll)[j], sep = '_'))
		file_j <- paste0(scratch_path, '/bbox_fiadhi_', layers[[ii]], '_', row.names(bbsll)[j], '.tif')
		point_stat[[ii]] <- statsByRadius(file_j, radii = radii, is_brick = FALSE)
		point_stat[[ii]]$variable <- layers[[ii]]
		if (file.exists(file_j)) deleteBox(file_j)
	}
	stats_by_point[[length(stats_by_point) + 1]] <- do.call('rbind', point_stat)
}	

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/geostats/dhi_',slice,'.r'))
