# Do arithmetic on rasters (sin, cos, and abs) and write to new rasters, within R.

library(raster)
sind <- function(n) sin(n * pi/180)
cosd <- function(n) cos(n * pi/180)

path_dem <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs'
path_tpi <- '/mnt/research/nasabio/data/tri_tpi'

aspect_big_in <- file.path(path_dem, 'conus_30m_aspect_big.tif')
aspect_small_in <- file.path(path_dem, 'conus_5k_aspect.tif')
tpi_big_in <- file.path(path_tpi, 'conus_30m_TPI_big.tif')
tpi_small_in <- file.path(path_tpi, 'conus_5k_dem_TPI.tif')

aspect_big <- raster(aspect_big_in)
aspect_big_sin <- sind(aspect_big)
aspect_big_cos <- cosd(aspect_big)
aspect_small <- raster(aspect_small_in)
aspect_small_sin <- sind(aspect_small)
aspect_small_cos <- cosd(aspect_small)

tpi_big <- raster(tpi_big_in)
tpi_big_abs <- abs(tpi_big)
tpi_small <- raster(tpi_small_in)
tpi_small_abs <- abs(tpi_small)

writeRaster(aspect_big_sin, file.path(path_dem, 'conus_30m_aspect_sin.tif'))
writeRaster(aspect_big_cos, file.path(path_dem, 'conus_30m_aspect_cos.tif'))
writeRaster(aspect_small_sin, file.path(path_dem, 'conus_5k_aspect_sin.tif'))
writeRaster(aspect_small_cos, file.path(path_dem, 'conus_5k_aspect_cos.tif'))

writeRaster(tpi_big_abs, file.path(path_tpi, 'conus_30m_TPI_abs.tif'))
writeRaster(tpi_small_abs, file.path(path_tpi, 'conus_5k_TPI_abs.tif'))