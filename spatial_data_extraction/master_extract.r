# GEODIVERSITY EXTRACTION BATCH SCRIPT v2
# Uses only GDAL commands for less memory usage.

# Edited 08 Jan 2018: Use $SCRATCH and $TMPDIR
# Edited 09 Jan 2018: Read directly from $SCRATCH, write to $TMPDIR
# Edited 11 Jan 2018: Get all variables from table. (rename to master_extract.r)
# Edited 09 Feb 2018: Update scratch path.
# Edited 19 Feb 2018: Add escape characters to file path.
# Edited 17 Apr 2018: Change BBS coordinates to the route midpoints, not centroids, and add smaller radii.

# Workflow:
# 1. Use extractBox() to make square of maximum radius (300 km) around focal point
# 2. Use extractFromCircle() to make a bunch of circles of different radii, getting summary statistics each time.

slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
tmp_path <- Sys.getenv('TMPDIR')
scratch_path <- '/mnt/ls15/scratch/groups/nasabio/VRTs'

# Boilerplate code to get the arguments passed in
args=(commandArgs(TRUE))

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

# Load functions to do the extraction
source('/mnt/research/nasabio/code/extractbox.r')
source('/mnt/research/nasabio/code/extractfromcircle_stack.r')

# Load variable table
vartable <- read.csv('/mnt/research/nasabio/data/geodiv_table_for_gdal.csv', stringsAsFactors = FALSE)

# Get all options from variable table.
k <- which(vartable$variable.id == geovar)
n_layers <- vartable$N.layers[k]
raster_file_name <- vartable$File.name[k]
max_radius <- vartable$Max.radius[k]
categorical <- geovar %in% c('gea','gea_5k','soil')

# Load coordinates for correct taxon.	
if (taxon == 'bbs') {
	n_slices <- vartable$N.slices.bbs[k]
	coords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv')
	radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500)
}	
if (taxon == 'fia') {
	n_slices <- vartable$N.slices.fiaall[k]
	source('/mnt/research/nasabio/code/loadfiaall.r')
	coords <- fiacoords
	radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500)
}

radii <- radii[radii <= max_radius]

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(coords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

stats_by_point <- list()

for (i in rowidxmin:rowidxmax) {
	print(i)
	
	focalpoint <- with(coords, cbind(lon, lat))[i,,drop=FALSE]	
	
	extractBox(coords = focalpoint,
			   raster_file = file.path(scratch_path, raster_file_name),
			   radius = max_radius,
			   fp = tmp_path,
			   filetags = paste(geovar, i, sep = '_'))

	file_i <- paste("\'", tmp_path, "\'", "/bbox_", geovar, "_", i, ".tif", sep="")
	
	stats_by_point[[length(stats_by_point) + 1]] <- 
		extractFromCircle(coords = focalpoint,
						  raster_file = file_i,
						  radii = radii,
						  fp = tmp_path,
						  filetag = paste(geovar, i, sep = '_'),
						  nlayers = n_layers,
						  is_categorical = categorical,
						  delete_temp = TRUE)
						  
	if (file.exists(paste(tmp_path, "/bbox_", geovar, "_", i, ".tif", sep=""))) {
		deleteBox(file_i)
	}
}

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/', taxon, '/allgeodiv_v2/', geovar, '_', slice, '.r'))
