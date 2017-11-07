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
  
#wcdem_albers <- projectRaster(wcdem, crs = aea_crs)
#writeRaster(wcdem_albers, file = '/mnt/research/nasabio/data/dem/SRTM_90m_DEM/west_coast_dem_albers.grd', format='raster')			  
			  
radii <- c(1000,5000,7500,10000,20000)

mns <- matrix(NA, nrow=nrow(fiaalbers), ncol=length(radii))
cvs <- mns			  

pb <- txtProgressBar(0, length(radii) * nrow(fiaalbers), style=3)
			  
for (r in 1:length(radii)) {
	for (i in 1:nrow(fiaalbers)) {
		setTxtProgressBar(pb, nrow(fiaalbers)*(r-1) + i)
		buf1 <- raster::buffer(fiaalbers[i,], width=r)
		buf1l <- spTransform(buf1, CRS=wgs_crs)
		# Error catcher:
		iserr <- try(extract(wcdem, fiaspatial[i,]), TRUE)
		if (!inherits(iserr, 'try-error')) {
			elevs <- na.omit(extract(wcdem, buf1l)[[1]])
			if (min(elevs) < 0) elevs <- elevs - min(elevs)
			mns[i,r] <- mean(elevs)
			cvs[i,r] <- sd(elevs)/mns[i,r]
		}
	}
}

close(pb)

# Convert to longform
elevdatalong <- data.frame(fiaspatial@coords, fiaspatial@data)
df2 <- data.frame(radius=rep(radii, each=nrow(elevdatalong)), mean_elev = as.vector(mns), cv_elev = as.vector(cvs))
elevdatalong <- cbind(elevdatalong, df2)
write.csv(elevdatalong, file='/mnt/research/nasabio/data/dem/SRTM_90m_DEM/elevstats90m.csv', row.names=FALSE)