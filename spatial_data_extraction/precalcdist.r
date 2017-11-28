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


###############################################

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


##############################################		
# Added 11 Sep: Precalculate distances for **all** combinations of radius and raster tile size.
# Use the same point for each.
# Edited 15 Sep: add NightLights and DHI rasters.
# Edited 28 Nov: add rasterToPoints() following suggestion on Stackoverflow.

source('/mnt/research/nasabio/code/extractbox.r')

makeDistRaster <- function(infile, outfile, radii, lon, lat) {
	require(raster)
	x <- raster(infile)
	
	# Get proximity map internally
	xpts <- rasterToPoints(x)
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

lon1 <- -98.59595
lat1 <- 39.79797

bbs_radii <- c(50, 75, 100, 150, 200, 300, 400, 500)
fia_radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500)

elevfile <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt'
bio1file <- '/mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.vrt'
bio5file <- '/mnt/research/nasabio/data/bioclim/Bioclim5k/rasterstack/bioclim5k_20170613.vrt'
geafile <- '/mnt/research/nasabio/data/geology/geo_ages/GEA.vrt'
stgfile <- '/mnt/research/nasabio/data/geology/soils/stg.vrt'
hffile <- '/mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf.vrt'
dhifparfile <- '/mnt/research/nasabio/data/dhi/dhi_fpar.vrt'
dhigppfile <- '/mnt/research/nasabio/data/dhi/dhi_gpp.vrt'
dhilai8file <- '/mnt/research/nasabio/data/dhi/dhi_lai8.vrt'
dhindvifile <- '/mnt/research/nasabio/data/dhi/dhi_ndvi.vrt'
nightlightfile <- '/mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl.vrt'

mat_names <- c('bio1','bio5','gea','stg','hf','dhifpar','dhigpp','dhilai8','dhindvi','nightlight')

for (n in mat_names) {
	print(n)
	nfile <- get(paste0(n, 'file'))
	print('Extracting box . . .')
	extractBox(coords = cbind(lon1, lat1),
			   raster_file = nfile,
			   radius = 500,
			   fp = '/mnt/research/nasabio/temp',
			   filetags = paste0(n ,'500'))
	print('Creating BBS matrix . . .')
	makeDistRaster(infile = paste0('/mnt/research/nasabio/temp/bbox_',n,'500.tif'),
				   radii = bbs_radii,
				   outfile = paste0('/mnt/research/nasabio/data/precalcdist/distlogical_bbs_',n,'.r'),
				   lon = lon1, lat = lat1)
	print('Creating FIA matrix . . .')
	makeDistRaster(infile = paste0('/mnt/research/nasabio/temp/bbox_',n,'500.tif'),
				   radii = fia_radii,
				   outfile = paste0('/mnt/research/nasabio/data/precalcdist/distlogical_fia_',n,'.r'),
				   lon = lon1, lat = lat1)
}
				   
extractBox(coords = cbind(lon1, lat1),
		   raster_file = elevfile,
		   radius = 300,
		   fp = '/mnt/research/nasabio/temp',
		   filetags = 'elev300')
makeDistRaster(infile = '/mnt/research/nasabio/temp/bbox_elev300.tif',
			   radii = bbs_radii[bbs_radii <= 300],
			   outfile = '/mnt/research/nasabio/data/precalcdist/distlogical_bbs_elev.r',
			   lon = lon1, lat = lat1)
makeDistRaster(infile = '/mnt/research/nasabio/temp/bbox_elev300.tif',
			   radii = fia_radii[fia_radii <= 300],
			   outfile = '/mnt/research/nasabio/data/precalcdist/distlogical_fia_elev.r',
			   lon = lon1, lat = lat1)		   