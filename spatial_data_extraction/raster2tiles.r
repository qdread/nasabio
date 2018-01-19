######################
# SPLIT RASTER
# See https://gis.stackexchange.com/questions/14712/splitting-raster-into-smaller-chunks-using-gdal

n_x <- 5
n_y <- 4

xs <- c(-125.0001, -64.99986)
ys <- c(19.99986, 56.00014)

xbounds <- seq(xs[1], xs[2], length.out = n_x + 1)
ybounds <- seq(ys[1], ys[2], length.out = n_y + 1)

input_file_name <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_aspect_big.vrt'

for (i in 1:n_x) {
	for (j in 1:n_y) {
		tile_dim <- paste('-projwin', xbounds[i], ybounds[j+1], xbounds[i+1], ybounds[j], '-a_srs "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"')
		output_file_name <- paste0('/mnt/research/nasabio/data/dem/SRTM_30m_DEM/aspect_tiles/aspect_', i, '_', j, '.tif')
		system2('gdal_translate', args = paste('-of GTIFF', tile_dim, input_file_name, output_file_name))
	}
}
