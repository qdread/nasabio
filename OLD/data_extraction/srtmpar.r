# Combine SRTM DEM into a single raster 
# Parallel version without doparallel.

n_slices <- 50
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))


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
  
# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fiaalbers),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]
			  
radii <- c(1000,5000,7500,10000,20000)

elevstats <- list()
			  
for (r in 1:length(radii)) {
	print(r)
	stats_r <- list()
	for (i in rowidxmin:rowidxmax) {
		buf1 <- raster::buffer(fiaalbers[i,], width=radii[r])
		buf1l <- spTransform(buf1, CRS=wgs_crs)
		# Error catcher:
		iserr <- try(extract(wcdem, buf1l), TRUE)
		if (!inherits(iserr, 'try-error')) {
			elevs <- na.omit(as.vector(extract(wcdem, buf1l)[[1]]))
			stats_r[[(i - rowidxmin + 1)]] <- c(mean_elev = mean(elevs), sd_elev = sd(elevs), min_elev=min(elevs), max_elev=max(elevs))
		}
		else {
		stats_r[[(i - rowidxmin + 1)]] <- c(mean_elev = NA, sd_elev = NA, min_elev=NA, max_elev=NA)
		}
	}
	elevstats[[length(elevstats) + 1]] <- do.call('rbind', stats_r)
}

elevstats <- do.call('rbind', elevstats)


# Convert to longform
elevdatalong <- data.frame(fiaspatial@coords[rowidxmin:rowidxmax, ], fiaspatial@data[rowidxmin:rowidxmax, ])
df2 <- data.frame(radius=rep(radii, each=nrow(elevdatalong)), mean_elev = elevstats[,1], sd_elev = elevstats[,2], min_elev = elevstats[,3], max_elev = elevstats[,4])
elevdatalong <- cbind(elevdatalong, df2)
write.csv(elevdatalong, file=paste0('/mnt/research/nasabio/data/fia/elevstats/elevstats90m_',slice,'.csv'), row.names=FALSE)