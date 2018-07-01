# GEODIVERSITY EXTRACTION BATCH SCRIPT v2: 1 kilometer edition
# QDR / NASAbioxgeo / 28 Jun 2018
# New version created for 1 kilometer radius only!
# Still has to be done in parallel unfortunately

# Workflow:
# 1. Use extractBox() to make square of maximum radius (300 km) around focal point
# 2. Use extractFromCircle() to make a bunch of circles of different radii, getting summary statistics each time.

slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
tmp_path <- Sys.getenv('TMPDIR')
scratch_path <- '/mnt/ffs17/groups/nasabio/VRTs'

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
	n_slices <- 1
	coords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv')
	radii <- c(1)
}	
if (taxon == 'fia') {
	n_slices <- 50
	source('/mnt/research/nasabio/code/loadfiaall.r')
	coords <- fiacoords
	radii <- c(1)
}

max_radius <- 1

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

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/', taxon, '/allgeodiv_1km/', geovar, '_', slice, '.r'))
