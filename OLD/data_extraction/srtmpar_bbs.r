# Extract elevation stats for BBS (t'other was fia)
# Parallel version without doparallel.

radii <- c(50000, 75000, 100000, 150000, 200000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
n_slices <- 50

combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]


# Extract values within a certain radius 
# See: https://rossijeanpierre.wordpress.com/2014/01/31/extracting-circular-buffers-from-a-raster-in-r-12/

library(raster)
library(dplyr)
library(rgeos)
library(rgdal)	

# Load Albers projection BBS coordinates (by route)

bbsalbers <- read.csv('/mnt/research/nasabio/data/bbs/bbs_aea_coords_byroute.csv')

# Fix NA problem.
bbsalbers[is.na(bbsalbers)] <- 0 # replace all NAs with 0,0 coordinate
		  
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

bbsalbers_sp <- SpatialPoints(coords = data.frame(x=bbsalbers$x_aea, y=bbsalbers$y_aea), proj4string = CRS(aea_crs))

# srtm 90m dem for west coast
wcdem <- raster('/mnt/research/nasabio/data/dem/SRTM_90m_DEM/west_coast_dem.vrt')
wcext <- extent(wcdem) 
 
# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(bbsalbers_sp@coords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]
			  
elevstats <- list()
			  
	for (i in rowidxmin:rowidxmax) {
		#print(i)
		buf1 <- raster::buffer(bbsalbers_sp[i,], width=r)
		buf1l <- spTransform(buf1, CRS=wgs_crs)
		# Check whether the buffer is inside the west coast
		bufext <- extent(buf1l)
		
		if (bufext@xmin > wcext@xmin & bufext@xmax < wcext@xmax & bufext@ymin > wcext@ymin & bufext@ymax < wcext@ymax) {
		
			# Error catcher:
			iserr <- try(extract(wcdem, buf1l), TRUE)
			if (!inherits(iserr, 'try-error')) {
				elevs <- na.omit(as.vector(extract(wcdem, buf1l)[[1]]))
				elevstats[[(i - rowidxmin + 1)]] <- c(mean_elev = mean(elevs), sd_elev = sd(elevs), min_elev=min(elevs), max_elev=max(elevs))
			}
			else {
				elevstats[[(i - rowidxmin + 1)]] <- c(mean_elev = NA, sd_elev = NA, min_elev=NA, max_elev=NA)
			}
		}
		else {
			elevstats[[(i - rowidxmin + 1)]] <- c(mean_elev = NA, sd_elev = NA, min_elev=NA, max_elev=NA)
		}
	}

elevstats <- do.call('rbind', elevstats)


# Convert to longform
elevdatalong <- data.frame(bbsalbers_sp@coords[rowidxmin:rowidxmax, ])
df2 <- data.frame(radius=r, mean_elev = elevstats[,1], sd_elev = elevstats[,2], min_elev = elevstats[,3], max_elev = elevstats[,4])
elevdatalong <- cbind(elevdatalong, df2)
write.csv(elevdatalong, file=paste0('/mnt/research/nasabio/data/bbs/elevstats/elevstats90m_', as.character(as.integer(r)),'_',slice,'.csv'), row.names=FALSE)