# Combine SRTM DEM into a single raster 

# Extract values within a certain radius 
# See: https://rossijeanpierre.wordpress.com/2014/01/31/extracting-circular-buffers-from-a-raster-in-r-12/

library(raster)
library(dplyr)
library(rgeos)
library(rgdal)	
		  
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)

fiacoords <- fiapnw %>%
	group_by(PLT_CN, STATECD, COUNTYCD, PLOT) %>%
	summarize(lat = mean(unique(LAT_FUZZSWAP)),
			  lon = mean(unique(LON_FUZZSWAP)))

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fiaspatial <- SpatialPointsDataFrame(coords = data.frame(x=fiacoords$lon, y=fiacoords$lat),
                                     data = as.data.frame(fiacoords[,1:4]),
                                     proj4string = CRS(wgs_crs)
)

fiaalbers <- spTransform(fiaspatial, CRSobj = CRS(aea_crs))		

# srtm 90m dem for west coast
wcdem <- raster('/mnt/research/nasabio/data/dem/SRTM_90m_DEM/west_coast_dem.vrt')
  
wcdem_albers <- projectRaster(wcdem, crs = aea_crs)
writeRaster(wcdem_albers, file = '/mnt/research/nasabio/data/dem/SRTM_90m_DEM/west_coast_dem_albers.grd', format='raster')			  
			  

# Test 5km buffer on one coordinate.

pbuf <- gBuffer(fiaalbers[100,], width=1000)

extract(wcdem_albers, pbuf)
