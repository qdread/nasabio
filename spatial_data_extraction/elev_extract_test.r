# Workflow:
# 1. Use extractBox() to make square of maximum radius (300 km) around focal point
# 1a. If needed, build VRT from this.
# 2. Use extractFromCircle() to make a bunch of circles of different radii, getting summary statistics each time.

# Test with one point

raster_file_name <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt'

focalpoint <- matrix(c(-100, 35), nrow = 1)
max_radius <- 300
scratch_path <- '/mnt/research/nasabio/temp'
radii <- c(50, 75, 100, 150, 200, 300)

source('/mnt/research/nasabio/code/extractbox.r')
source('/mnt/research/nasabio/code/extractfromcircle_stack.r')

extractBox(coords = focalpoint,
		   raster_file = raster_file_name,
		   radius = max_radius,
		   fp = scratch_path,
		   filetags = 'test')

extractFromCircle(coords = focalpoint,
				  raster_file = file.path(scratch_path, 'bbox_test.tif'),
				  radii = radii,
				  fp = scratch_path,
				  filetag = 'test',
				  nlayers = 1)

deleteBox(file.path(scratch_path, 'bbox_test.tif'))