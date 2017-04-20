
library(raster)
library(rgdal)	

usaelev <- raster('/mnt/research/nasabio/data/dem/usaelev.grd')
		  
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

usaelev_albers <- projectRaster(usaelev, crs = aea_crs)
writeRaster(usaelev_albers, file = '/mnt/research/nasabio/data/dem/usaelev_albers.grd', format='raster')			  