# Combine SRTM DEM into a single raster 

# Load rasters
setwd('/mnt/research/plz-lab/NEON/external_data/raw_external_data/srtm/')

x <- list()
for (i in dir(pattern = '*.tif')) x <- c(x, raster(i)) 

# Get raster extents
xextents <- lapply(x, extent)
xextents <- data.frame(xmin=sapply(xextents, '[', 1), xmax=sapply(xextents, '[', 2), ymin=sapply(xextents, '[', 3), ymax=sapply(xextents, '[', 4))

# 19 contains Puerto Rico, we can ignore it.
usaelev <- do.call(merge, x[1:18])

# Save the big raster.
writeRaster(usaelev, file = '/mnt/research/nasabio/data/dem/usaelev.grd', format='raster')

# Extract values within a certain radius 
# See: https://rossijeanpierre.wordpress.com/2014/01/31/extracting-circular-buffers-from-a-raster-in-r-12/

fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)

library(dplyr)
fiacoords <- fiapnw %>%
	group_by(PLT_CN, STATECD, COUNTYCD, PLOT) %>%
	summarize(lat = mean(unique(LAT_FUZZSWAP)),
			  lon = mean(unique(LON_FUZZSWAP)))

library(rgdal)			  
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fiaspatial <- SpatialPointsDataFrame(coords = data.frame(x=fiacoords$lon, y=fiacoords$lat),
                                     data = as.data.frame(fiacoords[,1:4]),
                                     proj4string = CRS(wgs_crs)
)

fiaalbers <- spTransform(fiaspatial, CRSobj = CRS(aea_crs))		  
usaelev_albers <- projectRaster(usaelev, crs = aea_crs)
writeRaster(usaelev_albers, file = '/mnt/research/nasabio/data/dem/usaelev_albers.grd', format='raster')			  
			  
library(rgeos)

# Test 5km buffer on one coordinate.
pbuf <- gBuffer(fiaalbers[1,], width=5000)
buf <- mask(usaelev_albers, pbuf)
buffer <- trim(buf, pad=2)