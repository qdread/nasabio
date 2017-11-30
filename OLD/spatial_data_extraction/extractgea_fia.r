
n_slices <- 1000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
raster_file_name <- '/mnt/research/nasabio/data/geology/geo_ages/GEA.vrt'
scratch_path <- Sys.getenv('TMPDIR')

# FIA lat long coordinates
source('/mnt/research/nasabio/code/loadfia.r')

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')

load('/mnt/research/nasabio/data/precalcdist/distlogical_fia_gea.r')

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
		   filetags = paste('fiagea', row.names(fiacoords)[j], sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_fiagea_', row.names(fiacoords)[j], '.tif')
	stats_by_point[[length(stats_by_point) + 1]] <- diversityByRadius(file_j, radii = radii, is_brick = FALSE)
	if (file.exists(file_j)) deleteBox(file_j)
}	

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/geostats/geoage_',slice,'.r'))
