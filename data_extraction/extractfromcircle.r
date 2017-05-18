# Function for extracting pixels from buffers around a focal point within a given radius, and calculating summary statistics
# This version uses Andy's code that does the following:
# make a polygon that approximates a circle and write it to a .shp
# make a system call to GDAL to clip the raster to the circle, writing it to a square .tif with missing values for pixels with centers outside the circle.

# issues to still resolve: 
# what to do about the fact that the raster is in lat long and the circle has a fixed radius in km? can we make a "great" circle?
# also, what about missing data?


# Author: QDR
# Project: NASABIOXGEO
# Date: 17 May 2017

# Needed inputs:
# matrix containing focal point coordinates (col1 is long, col2 is lat)
# path to raster containing the data (usually this will be the continental USA)
# radii is a vector of all the radii, in km, that we want to get summary stats for.
# fp is the path where the temporary files will be written



extract_from_buffer <- function(coords, raster_file, radii, lonbds, latbds, fp) {
	# define function to make an approximate circle with center point cx,cy, radius r, and n vertices
	approxCircle <- function(cx, cy, r, n){
		a <-  seq(0,  (2*pi), length.out=n+1)[-1]
		x <- cx + r * cos(a)
		y <- cy + r * sin(a)
		cbind(x, y)
	}
	
	# define function to make the polygon into a json
	polygon2json <- function(p){
		paste('{"type": "Polygon", "coordinates": [ [', paste(apply(p, 1, function(x){paste("[",x[1],", ",x[2],"]",sep="")}), collapse=","), '] ]}', sep="")
	}
	
	# Load raster and get its extent (bounding box).
	require(raster)
	ras <- raster(raster_file)
	ras_extent <- extent(ras)
	
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
			if (loni > ras_extent@xmin & loni < ras_extent@xmax & lati > ras_extent@ymin & lati < ras_extent@ymax) {
				
				for (r in 1:length(radii)) {
					# Create circle with 1000 vertices and create geojson object from it.
					circle_r <- approxCircle(loni, lati, radii[r], 1000)
					circle_json <- polygon2json(circle_r)
					circle_geojson <- readOGR(circle_json, "OGRGeoJSON", verbose = F,  p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
					# write the geojson to .shp
					writeOGR(circle_geojson, fp, "temp_circle", driver="ESRI Shapefile")
					# define arguments for system call
					call_args <- paste("-crop_to_cutline -overwrite -dstnodata NULL -cutline temp_circle.shp", raster_file, "temp_extracted.tif")
					# call GDAL twice, first to clip circle, then to calculate summary statistics on it.
					system2(command="gdalwarp", args=call_args)
					g.info <- system2(command="gdalinfo", args="-stats out.tif", stdout=TRUE)
					a <- g.info[(length(g.info)-3):length(g.info)]##get the last four lines
					summary_stats <- as.numeric(sapply(strsplit(a, split="="), "[[", 2)) # In order, this outputs: min, max, mean, sd
					stats_by_point[[i]][r, 2:5] <- summary_stats[c(3,4,1,2)]
				}
		
			}
		}
	}
	
	close(pb)
	return(stats_by_point)
	
}