# Function for extracting pixels from buffers around a focal point within a given radius, and calculating summary statistics
# gets the maximally sized square box, gets distances from each cell, then subsets by radius for any radius we want.
# I will also add subsampling, and possibly a way to account for partial pixels.

# Author: QDR
# Project: NASABIOXGEO
# Date: 15 May 2017

# Needed inputs:
# matrix containing focal point coordinates (col1 is long, col2 is lat)
# path to raster containing the data (usually this will be the continental USA)
# radii is a vector of all the radii, in km, that we want to get summary stats for.

extract_from_buffer <- function(coords, raster_file, radii) {
	
	# Load raster and get its extent (bounding box).
	require(raster)
	ras <- raster(raster_file)
	ras_extent <- extent(ras)
	
	# Initialize data structure to hold the summary stats for each radius at each point.
	# It is a list of data frames.
	# For now, it's just mean, sd, min, and max. More can be added if needed.
	emptydf <- data.frame(radius = radii, means = NA, sds = NA, mins = NA, maxes = NA, n = NA)
	stats_by_point <- replicate(nrow(coords), emptydf, simplify = FALSE)
	
	pb <- txtProgressBar(1, nrow(coords), style=3) # Tracks progress
	
	for (i in 1:nrow(coords)) {
		setTxtProgressBar(pb, i)
		
		# For each point, check whether the value is not missing and if the point is in the bounding box.
		if (!is.na(coords[i,1])) {
			loni <- coords[i,1]
			lati <- coords[i,2]
			if (loni > ras_extent@xmin & loni < ras_extent@xmax & lati > ras_extent@ymin & lati < ras_extent@ymax) {
				ei <- extract(ras, cbind(x=loni,y=lati), buffer=max(radii) * 1000, cellnumbers = TRUE) # The actual extraction of pixels, will be slow (buffer is in m)
				ex <- xFromCell(ras, ei[[1]][,1]) # Find longitudes of the extracted pixels (midpoints)
				ey <- yFromCell(ras, ei[[1]][,1]) # Find latitudes of the extracted pixels (midpoints)
				disti <- spDistsN1(cbind(ex,ey), c(loni,lati), longlat=TRUE) # Distance from each cell to point i
				
				for (r in 1:length(radii)) {
					eir <- na.omit(ei[[1]][disti <= radii[r], 2]) # Elevations within the circle
					stats_by_point[[i]][r, 2:6] <- c(mean(eir), sd(eir), min(eir), max(eir), length(eir)) # Calculate summary stats
				}
			}
		}
	}
	
	close(pb)
	return(stats_by_point)
	
}