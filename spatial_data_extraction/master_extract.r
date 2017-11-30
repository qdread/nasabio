# MASTER GEODIVERSITY EXTRACTION SCRIPT
# Works with any taxon and geodiversity variable
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
scratch_path <- Sys.getenv('TMPDIR')

# Boilerplate code to get the arguments passed in
args=(commandArgs(TRUE))

for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')	
	
# Load variable table
vartable <- read.csv('/mnt/research/nasabio/data/geodiv_stat_table_separate.csv', stringsAsFactors = FALSE)

# Get all options from variable table.
k <- which(vartable$variable.id == geovar)
n_layers <- vartable$N.layers[k]
raster_file_name <- file.path('/mnt/research/nasabio/data', vartable$File.location[k], vartable$File.name[k])
max_radius <- vartable$Max.radius[k]
use_abs <- vartable$Use.abs.value[k] == 'yes'
use_trig <- vartable$Use.trig[k] == 'yes'
use_categorical <- vartable$Use.categorical[k] == 'yes'	
	
# Load coordinates for correct taxon.	
if (taxon == 'bbs') {
	n_slices <- vartable$N.slices.bbs[k]
	coords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')
	radii <- c(50, 75, 100, 150, 200, 300, 400, 500)
}	
if (taxon == 'fia') {
	n_slices <- vartable$N.slices.fiapnw[k]
	source('/mnt/research/nasabio/code/loadfia.r')
	coords <- fiacoords
	radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500)
}

radii <- radii[radii <= max_radius]

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(coords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Call extraction function, specifying the raster from which to extract data.
stats_by_point <- list()

if (rowidxmin < rowidxmax) {
	pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)
}
		   
for (j in rowidxmin:rowidxmax) {
	if (rowidxmin < rowidxmax) { 
		setTxtProgressBar(pb, j) 
	}
	focalpoint <- with(coords, cbind(lon, lat))[j,,drop=FALSE]	
	extractBox(coords = focalpoint,
			   raster_file = raster_file_name,
			   radius = max_radius,
			   fp = scratch_path,
			   filetags = paste(taxon, geovar, row.names(coords)[j], sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_', taxon, '_', geovar, '_', row.names(coords)[j], '.tif')
	stats_by_point[[length(stats_by_point) + 1]] <- 
		statsByRadius(boxfile = file_j, 
					  centerpoint = focalpoint, 
					  radii = radii, 
					  is_brick = n_layers > 1, 
					  is_abs = use_abs, 
					  is_trig = use_trig, 
					  is_categorical = use_categorical)
	if (file.exists(file_j)) {
		deleteBox(file_j)
	}
}	

if (rowidxmin < rowidxmax) { 
	close(pb) 
}

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/', taxon, '/allgeodiv/', geovar, '_', slice, '.r'))
