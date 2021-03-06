library(raster)
library(rgeos)
library(rgdal)	
wcdem <- raster('/mnt/research/nasabio/data/dem/SRTM_90m_DEM/west_coast_dem.vrt')
wcdem_albers <- projectRaster(wcdem, crs = aea_crs)
writeRaster(wcdem_albers, file = '/mnt/research/nasabio/data/dem/SRTM_90m_DEM/west_coast_dem_albers.grd', format='raster')			  
