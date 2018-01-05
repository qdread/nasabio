# ELEVATION EXTRACTION BATCH SCRIPT
# Uses only GDAL commands for less memory usage.

# Workflow:
# 1. Use extractBox() to make square of maximum radius (300 km) around focal point
# 2. Use extractFromCircle() to make a bunch of circles of different radii, getting summary statistics each time.

slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
scratch_path <- Sys.getenv('TMPDIR')
n_slices <- 10000

# Boilerplate code to get the arguments passed in
args=(commandArgs(TRUE))

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

raster_file_names <- c(elevation = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt',
					   slope = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_slope_big.vrt',
					   roughness = '/mnt/research/nasabio/data/tri_tpi/conus_30m_dem_roughness.vrt',
					   tpi = '/mnt/research/nasabio/data/tri_tpi/conus_30m_dem_TPI.vrt',
					   tri = '/mnt/research/nasabio/data/tri_tpi/conus_30m_dem_TRI.vrt')

raster_file_name <- raster_file_names[geovar]
max_radius <- 300
radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300)

source('/mnt/research/nasabio/code/extractbox.r')
source('/mnt/research/nasabio/code/extractfromcircle_stack.r')
source('/mnt/research/nasabio/code/loadfiaall.r')

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

stats_by_point <- list()

for (i in rowidxmin:rowidxmax) {
	print(i)
	
	focalpoint <- with(fiacoords, cbind(lon, lat))[i,,drop=FALSE]	
	
	extractBox(coords = focalpoint,
			   raster_file = raster_file_name,
			   radius = max_radius,
			   fp = scratch_path,
			   filetags = paste(geovar, i, sep = '_'))

	file_i <- file.path(scratch_path, paste0('bbox_', geovar, '_', i, '.tif'))		   
			   
	stats_by_point[[length(stats_by_point) + 1]] <- 
		extractFromCircle(coords = focalpoint,
						  raster_file = file_i,
						  radii = radii,
						  fp = scratch_path,
						  filetag = paste(geovar, i, sep = '_'),
						  nlayers = 1)
						  
	if (file.exists(file_i)) {
		deleteBox(file_i)
	}
}

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats_usa/', geovar, '_', slice, '.r'))
