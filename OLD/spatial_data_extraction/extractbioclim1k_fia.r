
n_slices <- 1000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
raster_file_name <- '/mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.vrt'
scratch_path <- Sys.getenv('TMPDIR')

# FIA lat long coordinates
source('/mnt/research/nasabio/code/loadfia.r')

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')

load('/mnt/research/nasabio/data/precalcdist/distlogical_fia_bio1.r')

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Radii of circles where we want summary stats
radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500) # in km

# Call extraction function, specifying the raster from which to extract data.
# 1 km bioclim
stats_by_point <- list()

pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)
		   
for (j in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, j)
	extractBox(coords = with(fiacoords, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = raster_file_name,
		   radius = 500,
		   fp = scratch_path,
		   filetags = paste('fia1k', row.names(fiacoords)[j], sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_fia1k_', row.names(fiacoords)[j], '.tif')
	stats_by_point[[length(stats_by_point) + 1]] <- statsByRadius(file_j, radii = radii, is_brick = TRUE)
	if (file.exists(file_j)) deleteBox(file_j)
}	

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/climstats/bioclim1k_',slice,'.r'))
