# Precalculate distances for each cell the large bounding box 300 km radius.
# This will save a lot of time. We can just extract the cell numbers required each time.

fp <- '/mnt/research/nasabio/temp'

library(raster)
x <- raster(file.path(fp, 'bbox_1_fiatest.tif'))

# Get proximity map internally
xvals <- as.data.frame(x)
xcoords <- xyFromCell(x, 1:ncell(x))
xdist <- spDistsN1(xcoords, c(loni,lati), longlat=TRUE)

radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300)

idx_r <- list()

for (r in 1:length(radii)) {
	idx_r[[r]] <- xdist <= radii[r]
}

# Create logical matrix.
# Use column r to subset xvalues each time.
idx_r <- do.call('cbind', idx_r)

save(idx_r, file = '/mnt/research/nasabio/data/fia/distlogical.r')

####################################################

# general function

makeDistRaster <- function(infile, outfile, radii, lon, lat) {
	require(raster)
	x <- raster(infile)
	
	# Get proximity map internally
	xvals <- as.data.frame(x)
	xcoords <- xyFromCell(x, 1:ncell(x))
	xdist <- spDistsN1(xcoords, c(lon,lat), longlat=TRUE)
	
	idx_r <- list()
	
	for (r in 1:length(radii)) {
		idx_r[[r]] <- xdist <= radii[r]
	}

	# Create logical matrix.
	# Use column r to subset xvalues each time.
	idx_r <- do.call('cbind', idx_r)

	save(idx_r, file = outfile)

}


###############################################3

# 1 km tile, 500 km radius.

bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')
raster_file_name <- '/mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.vrt'
source('/mnt/research/nasabio/code/extractbox.r')

j<-1
extractBox(coords = with(bbsll, cbind(lon, lat))[j,,drop=FALSE],
	   raster_file = raster_file_name,
	   radius = 500,
	   fp = '/mnt/research/nasabio/temp',
	   filetags = paste('bbstest', row.names(bbsll)[j], sep = '_'))
	   
makeDistRaster(infile = '/mnt/research/nasabio/temp/bbox_bbstest_1.tif',
			   radii = c(50, 75, 100, 150, 200, 300, 400, 500),
			   outfile = '/mnt/research/nasabio/data/fia/distlogical_1ktile.r',
			   lon = bbsll$lon[j], lat = bbsll$lat[j])

# 5 km tile, 500 km radius.
			   
extractBox(coords = with(bbsll, cbind(lon, lat))[j,,drop=FALSE],
	   raster_file = '/mnt/research/nasabio/data/bioclim/Bioclim5k/rasterstack/bioclim5k_20170613.vrt',
	   radius = 500,
	   fp = '/mnt/research/nasabio/temp',
	   filetags = paste('bbstest5k', row.names(bbsll)[j], sep = '_'))

makeDistRaster(infile = '/mnt/research/nasabio/temp/bbox_bbstest5k_1.tif',
			   radii = c(50, 75, 100, 150, 200, 300, 400, 500),
			   outfile = '/mnt/research/nasabio/data/fia/distlogical_5ktile.r',
			   lon = bbsll$lon[j], lat = bbsll$lat[j])
	   