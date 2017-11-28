# Extract box that always has the same number of cells regardless of what latitude.

extractBox <- function(lon, lat, radius, input_file, input_proj, output_proj, output_res, output_file_path, output_file_name, lonbds = c(-125, -67), latbds = c(25, 50)) {
	require(sp)
	require(rgdal)
	require(raster)
	
	if (is.na(lon) | is.na(lat)) {
		return("Missing coordinate.")
	}
	
	if (lon < lonbds[1] | lon > lonbds[2] | lat < latbds[1] | lat > latbds[2]) {
		return("Point is not inside area of interest.")
	}
	
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
	
	# Calculate lat and long boundaries of box
	box_r <- latLongBox(lon, lat, radius)
	
	# make box to spatial object
	b_box <- SpatialPoints(coords = box_r, proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
	# Project lat long of box into the raster's CRS
	b_box <- spTransform(b_box, CRSobj = CRS(input_proj))
	
	# Create geojson object from bounding box
	box_json <- polygon2json(as(extent(b_box), 'SpatialPolygons')@polygons[[1]]@Polygons[[1]]@coords)
	box_geojson <- readOGR(box_json, "OGRGeoJSON", verbose = FALSE,  p4s = input_proj)
	# write the geojson to .shp
	writeOGR(box_geojson, output_file_path, "temp_bbox", driver="ESRI Shapefile")
	
	# define arguments for system call
	output_proj <- paste0('"', output_proj, '"')
	call_args <- paste("-t_srs", output_proj, "-tr", output_res, output_res, "-crop_to_cutline -overwrite -dstnodata NULL -cutline", file.path(output_file_path, "temp_bbox.shp"), input_file, file.path(output_file_path,output_file_name))
	# call GDAL to clip the box.
	system2(command="gdalwarp", args=call_args)	
	
	# Remove temporary shapefile
	for (file_ext in c('.shp', '.shx', '.prj', '.dbf')) system2(command="rm", args=file.path(output_file_path, paste0("temp_bbox", file_ext)))
}
		