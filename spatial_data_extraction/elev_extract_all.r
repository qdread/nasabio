# ELEVATION EXTRACTION BATCH SCRIPT
# Uses only GDAL commands for less memory usage.

# Edited 08 Jan 2018: Use $SCRATCH and $TMPDIR

# Workflow:
# 1. Use extractBox() to make square of maximum radius (300 km) around focal point
# 2. Use extractFromCircle() to make a bunch of circles of different radii, getting summary statistics each time.

slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
tmp_path <- Sys.getenv('TMPDIR')
scratch_path <- Sys.getenv('SCRATCH')
n_slices <- 10000

# Boilerplate code to get the arguments passed in
args=(commandArgs(TRUE))

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

# Names of vrts
raster_file_names <- c(elevation = 'dem/conus_30m_dem_big_singlefile.vrt',
					   slope = 'dem/VRTs/conus_30m_slope_big.vrt',
					   roughness = 'tri_tpi/conus_30m_dem_roughness.vrt',
					   tpi = 'tri_tpi/conus_30m_dem_TPI.vrt',
					   tri = 'tri_tpi/conus_30m_dem_TRI.vrt'))

# Names of tifs
raw_file_names <- c(elevation = 'dem/conus_30m_dem_big.tif',
				    slope = 'dem/conus_30m_slope_big.tif',
				    roughness = 'tri_tpi/conus_30m_roughness_big.tif',
				    tpi = 'tri_tpi/conus_30m_TPI_big.tif',
				    tri = 'tri_tpi/conus_30m_TRI_big.tif'))				   
					   
# Copy VRT and its corresponding TIF to TMPDIR.
raster_file_name <- raster_file_names[geovar]
raw_file_name <- raw_file_names[geovar]
system2('cp', args = paste(file.path(scratch_path, raster_file_name), tmp_path)
system2('cp', args = paste(file.path(scratch_path, raw_file_name), tmp_path)
raster_file_name <- file_path(tmp_path, raster_file_name)

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
			   fp = tmp_path,
			   filetags = paste(geovar, i, sep = '_'))

	file_i <- file.path(tmp_path, paste0('bbox_', geovar, '_', i, '.tif'))		   
			   
	stats_by_point[[length(stats_by_point) + 1]] <- 
		extractFromCircle(coords = focalpoint,
						  raster_file = file_i,
						  radii = radii,
						  fp = tmp_path,
						  filetag = paste(geovar, i, sep = '_'),
						  nlayers = 1)
						  
	if (file.exists(file_i)) {
		deleteBox(file_i)
	}
}

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats_usa/', geovar, '_', slice, '.r'))
