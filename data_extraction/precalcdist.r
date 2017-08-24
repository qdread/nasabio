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