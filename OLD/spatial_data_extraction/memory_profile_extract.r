# Memory usage profile: rasters
# Create rasters of different sizes and resolutions, calculate basic functions on them.

# Source functions
# ----------------

source('/mnt/research/nasabio/code/extractbox.r')

# Memory profiling
# ----------------

library(pryr)
library(lineprof)

# Test the code with a random point in the USA, 300 km radius
focalpoint <- matrix(c(-100, 35), nrow = 1)
max_radius <- 300
scratch_path <- '/mnt/research/nasabio/temp'
radii <- c(50, 75, 100, 150, 200, 300, 400, 500)

# Track memory usage for extracting the box and loading pixel values into memory
profile_extract <- lineprof(
	extractBox(coords = focalpoint,
		   raster_file = raster_file_name,
		   radius = max_radius,
		   fp = scratch_path,
		   filetags = 'test')
)

file_name <- paste0(scratch_path, '/bbox_', 'test', '.tif')
 
# Track memory usage for calculating summary statistics on the pixels
profile_calc_stats <- lineprof(
	statsByRadius(boxfile = file_name, 
				  centerpoint = focalpoint, 
				  radii = radii, 
				  is_brick = FALSE, 
				  is_abs = FALSE, 
				  is_trig = FALSE, 
				  is_categorical = FALSE)
)

deleteBox(file_j)

