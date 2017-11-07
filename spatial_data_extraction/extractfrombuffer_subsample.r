# Function for extracting pixels from buffers around a focal point within a given radius, and calculating summary statistics
# gets the maximally sized square box, gets distances from each cell, then subsets by radius for any radius we want.
# Version including subsampling.

# Author: QDR
# Project: NASABIOXGEO
# Date: 15 May 2017

# Function for sampling n random points inside a circle of radius r.
# See http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
rcircle <- function(n=1, r=1) {
	a1 <- runif(n)
	b1 <- runif(n)
	a <- pmin(a1, b1)
	b <- pmax(a1, b1)
	cbind(b * r * cos(2 * pi * a / b), b * r * sin(2 * pi * a / b))
}

# Needed inputs:
# matrix containing focal point coordinates (col1 is long, col2 is lat)
# path to raster containing the data (usually this will be the continental USA)
# radius is the single radius, in km, that we want to get summary statistics for.
# samplesize is size, in number of pixels, of the subsample to be taken within the buffer.

extract_from_buffer_subsample <- function(coords, raster_file, radius, samplesize) {
	
	# Load raster and get its extent (bounding box).
	require(raster)
	ras <- raster(raster_file)
	ras_extent <- extent(ras)
	
	# Initialize data structure to hold the summary stats for each radius at each point.
	# It is a data frame.
	# For now, it's just mean, sd, min, and max. More can be added if needed.
	stats_by_point <- data.frame(lon = coords[,1], lat = coords[,2], radius = radius, samplesize = samplesize, mean = NA, sd = NA, min = NA, max = NA, n = 0)
	
	pb <- txtProgressBar(1, nrow(coords), style=3) # Tracks progress
	
	# proj4 strings for long lat and for our chosen projection.
	aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
	wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'
		
	for (i in 1:nrow(coords)) {
		setTxtProgressBar(pb, i)
		
		# For each point, check whether the value is not missing and if the point is in the bounding box.
		if (!is.na(coords[i,1])) {
			loni <- coords[i,1]
			lati <- coords[i,2]
			if (loni > ras_extent@xmin & loni < ras_extent@xmax & lati > ras_extent@ymin & lati < ras_extent@ymax) {
							
				# Convert the focal point coordinate to a projected coordinate system
				albersi <- spTransform(SpatialPoints(cbind(loni, lati), proj4string=CRS(wgs_crs)), CRSobj = CRS(aea_crs))
				
				# Sample samplesize points inside a circle of the given radius.
				randompts <- rcircle(n = samplesize, r = radius * 1000)
				# Add the coordinate origin and convert back to longlat.
				randompts <- randompts + t(replicate(samplesize, albersi@coords[1,]))
				extractpts <- spTransform(SpatialPoints(randompts, proj4string=CRS(aea_crs)), CRSobj = CRS(wgs_crs))
			
				ei <- extract(ras, extractpts) # The actual extraction of pixels, will be slow
				ei_nona <- na.omit(as.vector(ei))    # Remove missing values.
				if (length(ei_nona) > 0) stats_by_point[i, 5:9] <- c(mean(ei_nona), sd(ei_nona), min(ei_nona), max(ei_nona), length(ei_nona)) # Calculate summary stats
			}
		}
	}
	
	
	close(pb)
	return(stats_by_point)
	
}