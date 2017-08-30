# Function to loop through a coordinate list and a raster file
# At each iteration, extract a square with given radius and save to .hdf somewhere

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

# Use precalculated distances to get the stats for the pixels within each radius.
# idx_r should be loaded
statsByRadius <- function(boxfile) {
	radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300)
	if (!file.exists(boxfile)) return(NA)
	x <- raster(boxfile)
	xvals <- as.data.frame(x)
	stats_r <- list()
	for (r in 1:ncol(idx_r)) {
		vals_r <- xvals[idx_r[, r], , drop = FALSE]
		means <- apply(vals_r, 2, mean, na.rm = TRUE)
		sds <- apply(vals_r, 2, sd, na.rm = TRUE)
		mins <- apply(vals_r, 2, min, na.rm = TRUE)
		maxes <- apply(vals_r, 2, max, na.rm = TRUE)
		nums <- apply(vals_r, 2, function(z) sum(!is.na(z)))
		stats_r[[r]] <- data.frame(radius = radii[r],
								   variable = names(xvals),
								   mean = means,
								   sd = sds,
								   min = mins,
								   max = maxes,
								   n = nums)
	}
	do.call('rbind', stats_r)	
}

# delete pre-created tif
deleteBox <- function(boxfile) system2("rm", args=boxfile)
