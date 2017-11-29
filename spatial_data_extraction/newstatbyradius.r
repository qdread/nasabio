# New stats by radius workflow

# Must define radii (vector of n different radii) and focalpoint (vector of length 2: lon, lat)

x <- raster('/mnt/research/nasabio/temp/bbox_elev300.tif') # Load raster
xpts <- rasterToPoints(x) # Convert raster to matrix of points
xdists <- spDistsN1(pts = xpts[,1:2], pt = focalpoint, longlat = TRUE) # Calculate distance to center
xsubsets <- sapply(radii, function(z) xdists <= z) # Logicals for each radius circle
xdat <- apply(xsubsets, 2, function(z) xpts[z,3]) # Data indexed by the subsets
xsd <- lapply(xdat, sd, na.rm = TRUE) # Example summary statistic calculation

# Also try to do this with bricks so that we can extract multiple rasters, make them into a brick, then do the distances and stats all at once to save time.

namepart1 <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_'
namepart2 <- '_big.vrt'
varname <- c('dem','aspect','slope','TPI')

source('/mnt/research/nasabio/code/extractbox.r')

focalpoint <- c(-120, 37.5)

# Extract the four different topography stats
for (k in varname) {
extractBox(coords = matrix(focalpoint, nrow=1),
		   raster_file = paste0(namepart1, k, namepart2),
		   radius = 300,
		   fp = '/mnt/research/nasabio/temp',
		   filetags = paste('test', k, sep = '_'))
}

# Combine them into a stack (bricks only used for multiple files)

fnames <- paste0('/mnt/research/nasabio/temp/bbox_test_', varname, '.tif')

xstack <- stack(as.list(fnames))
system.time({xpts <- rasterToPoints(xstack)})
# This is way too slow. Should either make a brick or do them all separately.