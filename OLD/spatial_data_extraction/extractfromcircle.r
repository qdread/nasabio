# Function for extracting pixels from buffers around a focal point within a given radius, and calculating summary statistics
# Author: QDR
# Project: NASABIOXGEO
# Date: 17 May 2017

# Last modified: 04 Aug 2017 (added some config options to gdalwarp, attempting to make parallel read possible)
# Last modified: 22 Jun 2017 (changed tif to hdf)
# Last modified: 18 May 2017 (made circle into a "great" circle using the proper spherical geometry with lat-long coordinates.)

# Needed inputs:
# matrix containing focal point coordinates (col1 is long, col2 is lat)
# path to raster containing the data (usually this will be the continental USA)
# radii is a vector of all the radii, in km, that we want to get summary stats for.
# lonbds and latbds are each a vector of the min and max longitude and latitude. Default is continental usa
# fp is the path where the temporary files will be written
# filetag is a tag that prevents overwriting of files.

# What it does:
# Make a polygon approximating a circle and write it to a .shp
# make a system call to GDAL to clip the raster to the circle, writing it to a square .tif with missing values for pixels with centers outside the circle.
# call GDAL again to get summary stats and parse the GDAL output file in R.

extractFromCircle <- function(coords, raster_file, radii, lonbds = c(-125, -67), latbds = c(25, 50), fp, filetag = '') {
	require(sp)
	require(rgdal)

	# define function to make an approximate circle with center point lon1, lat1, radius d, and n vertices
	latLongCircle <- function(lon1, lat1, d, n = 1000) {
  
	  lat1 <- lat1 * pi/180
	  lon1 <- lon1 * pi/180
	  d <- d/6371			# Earth's radius in km
	  
	  thetas <- seq(0, 2*pi, length.out = n)
	  lat <- asin(sin(lat1)*cos(d)+cos(lat1)*sin(d)*cos(thetas))
	  dlon <- atan2(sin(thetas)*sin(d)*cos(lat1), cos(d)-sin(lat1)*sin(lat))
	  lon = ((lon1+dlon+pi) %% (2*pi)) - pi
	  
	  cbind(lon, lat) * 180/pi
  
	}
	
	# define function to make the polygon into a json
	polygon2json <- function(p){
		paste('{"type": "Polygon", "coordinates": [ [', paste(apply(p, 1, function(x){paste("[",x[1],", ",x[2],"]",sep="")}), collapse=","), '] ]}', sep="")
	}
	
	# Initialize data structure to hold the summary stats for each radius at each point.
	# It is a list of data frames.
	# For now, it's just mean, sd, min, and max. More can be added if needed.
	emptydf <- data.frame(radius = radii, mean = NA, sd = NA, min = NA, max = NA)
	stats_by_point <- replicate(nrow(coords), emptydf, simplify = FALSE)
	
	pb <- txtProgressBar(1, nrow(coords), style=3) # Tracks progress
	
	for (i in 1:nrow(coords)) {
		setTxtProgressBar(pb, i)
		
		
		# For each point, check whether the value is not missing and if the point is in the bounding box.
		if (!is.na(coords[i,1])) {
			loni <- coords[i,1]
			lati <- coords[i,2]
			if (loni > lonbds[1] & loni < lonbds[2] & lati > latbds[1] & lati < latbds[2]) {
				
				for (r in 1:length(radii)) {
					# Create circle with 1000 vertices and create geojson object from it.
					circle_r <- latLongCircle(loni, lati, radii[r], 1000)
					circle_json <- polygon2json(circle_r)
					circle_geojson <- readOGR(circle_json, "OGRGeoJSON", verbose = F,  p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
					# write the geojson to .shp
					writeOGR(circle_geojson, fp, paste0("temp_circle_", filetag), driver="ESRI Shapefile")
					# define arguments for system call
					call_args <- paste("--config VRT_SHARED_SOURCE 0 --config GDAL_MAX_DATASET_POOL_SIZE 1024 -crop_to_cutline -overwrite -dstnodata NULL -cutline", file.path(fp, paste0("temp_circle_", filetag, ".shp")), raster_file, file.path(fp,paste0("temp_circle_extracted", filetag, ".hdf")))
						
					# call GDAL twice, first to clip circle, then to calculate summary statistics on it.
					system2(command="gdalwarp", args=call_args)
					g.info <- system2(command="gdalinfo", args=paste("-stats", file.path(fp,paste0("temp_circle_extracted", filetag, ".hdf"))), stdout=TRUE) 
					
					a <- g.info[(length(g.info)-3):length(g.info)] # get the last four lines
					summary_stats <- as.numeric(sapply(strsplit(a, split="="), "[[", 2)) # In order, this outputs: min, max, mean, sd
					stats_by_point[[i]][r, 2:5] <- summary_stats[c(3,4,1,2)] # write summary stats to the output list.
					
					# Delete the temporarily created files between each iteration because GDAL gets mad if you overwrite an existing file.
					system2(command="rm", args=file.path(fp, paste0("temp_circle*", filetag, "*")))
				}
		
			}
		}
	}
	
	close(pb)
	return(stats_by_point)
	
}