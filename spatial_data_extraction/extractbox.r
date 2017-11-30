# Function to loop through a coordinate list and a raster file
# At each iteration, extract a square with given radius and save to .tif
# Edited 30 Nov: new stats function that does not use precalculated distance matrix.
# Edited 10 Nov: add option to get mean of absolute value (used for TPI)

extractBox <- function(coords, raster_file, radius, lonbds = c(-125, -67), latbds = c(25, 50), fp, filetags = 1:nrow(coords), progress = FALSE) {
	require(sp)
	require(rgdal)
	
	# define function to make box centered at a given coordinate with a given radius
	# using lat long
	
	latLongBox <- function(lon1, lat1, r) {
		# length of 1 degree latitude anywhere, and 1 degree longitude at equator, in km
		length1 <- 6371*2*pi/360
		latmin <- lat1 - r/length1
		latmax <- lat1 + r/length1
		lonmin <- lon1 - cos(lat1 * pi/180) * r/length1
		lonmax <- lon1 + cos(lat1 * pi/180) * r/length1
		cbind(c(lonmin, lonmin, lonmax, lonmax), c(latmin, latmax, latmax, latmin)) 
	}
	
	# define function to make the polygon into a json
	polygon2json <- function(p) {
		paste('{"type": "Polygon", "coordinates": [ [', paste(apply(p, 1, function(x){paste("[",x[1],", ",x[2],"]",sep="")}), collapse=","), '] ]}', sep="")
	}
	
	if (progress) pb <- txtProgressBar(1, nrow(coords), style=3) # Tracks progress
	
	for (i in 1:nrow(coords)) {
		if (progress) setTxtProgressBar(pb, i)
				
		# For each point, check whether the value is not missing and if the point is in the bounding box.
		if (!is.na(coords[i,1])) {
			loni <- coords[i,1]
			lati <- coords[i,2]
			if (loni > lonbds[1] & loni < lonbds[2] & lati > latbds[1] & lati < latbds[2]) {
				# Create bounding box and create geojson object from it.
				box_r <- latLongBox(loni, lati, radius)
				box_json <- polygon2json(box_r)
				box_geojson <- readOGR(box_json, "OGRGeoJSON", verbose = FALSE,  p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
				# write the geojson to .shp
				writeOGR(box_geojson, fp, paste0("temp_bbox_", filetags[i]), driver="ESRI Shapefile")
				# define arguments for system call
				call_args <- paste("-crop_to_cutline -overwrite -dstnodata NULL -cutline", file.path(fp, paste0("temp_bbox_", filetags[i], ".shp")), raster_file, file.path(fp,paste0("bbox_", filetags[i], ".tif")))
				# call GDAL to clip the box.
				system2(command="gdalwarp", args=call_args)
				# Remove temporary shapefile
				for (file_ext in c('.shp', '.shx', '.prj', '.dbf')) system2(command="rm", args=file.path(fp, paste0("temp_bbox_", filetags[i], file_ext)))
			}
		}
	}
	
	if (progress) close(pb)
}

# Definition of which summary stats we want to calculate
# Add more here if you want.
# Only summary stats that ignore the spatial arrangement of the points can be calculated here.
summstats_continuous <- function(z) {
	z <- na.omit(z)
	c(mean = mean(z),
	  sd = sd(z),
	  min = min(z),
	  max = max(z),
	  n = length(z))
}

summstats_categorical <- function(z) {
	z <- na.omit(z)
	ztable <- table(z)
	zprop <- ztable/sum(ztable)
	c(richness = length(unique(z)),
	  diversity = -sum(zprop * log(zprop)),
	  n = length(z))
}

# Edited stats function 30 Nov.
# Unfortunately the distance matrix must be calculated here.
# This function now works for both categorical and continuous data.

statsByRadius <- function(boxfile, 
						  centerpoint, 
						  radii, 
						  is_brick = FALSE, 
						  is_trig = FALSE, 
						  is_abs = FALSE, 
						  is_categorical = FALSE) {
	require(raster)
	if (!file.exists(boxfile)) {
		return(NA)
	}
	if (!is_brick) {
		x <- raster(boxfile)
	} else {
		x <- brick(boxfile)
	}
	xpts <- rasterToPoints(x) # Convert raster to matrix of points
	xdists <- spDistsN1(pts = xpts[, 1:2], pt = matrix(centerpoint, nrow=1), longlat = TRUE) # Calculate distance to center
	xsubsets <- sapply(radii, function(z) xdists <= z) # Logicals for each radius circle
	xvals <- apply(xsubsets, 2, function(z) xpts[z, -(1:2), drop = FALSE]) # Data indexed by the subsets
	
	# Calculation of summary statistics
	
	if (!is_categorical) {
		# continuous data
		if (!is_trig) {
			if (is_abs) {
				# convert to absolute value before calculating statistics
				xvals <- lapply(xvals, abs)
			}
			# calculate summary statistics
			xstats <- do.call('rbind', lapply(xvals, function(y) t(apply(y, 2, summstats_continuous))))
		} else {
			# find sin and cos and calculate stats on them (degrees to radians first)
			xsin <- lapply(xvals, function(z) sin(pi/180 * z))
			xcos <- lapply(xvals, function(z) cos(pi/180 * z))
			sinstats <- lapply(xsin, function(y) t(apply(y, 2, summstats_continuous)))
			cosstats <- lapply(xcos, function(y) t(apply(y, 2, summstats_continuous)))
			# combine sin and cos summary stats into single data frame
			xstats <- do.call('rbind', mapply(cbind, sinstats, cosstats, SIMPLIFY = FALSE))
			dimnames(xstats)[[2]] <- paste(dimnames(xstats)[[2]], rep(c('sin', 'cos'), each = length(dimnames(xstats)[[2]])/2), sep = '_')
		}
	} else {
		# categorical data
		xstats <- do.call('rbind', lapply(xvals, function(y) t(apply(y, 2, summstats_categorical))))
	}
	# Add radius and variable information to the output data
	data.frame(radius = rep(radii, each = ncol(xvals[[1]])),
			   variable = dimnames(xstats)[[1]],
			   xstats)
}

# delete pre-created tif
deleteBox <- function(boxfile) system2("rm", args=boxfile)
