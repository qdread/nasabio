# Crop Michigan from different rasters and save the results as new, smaller rasters.
# Edit: added arguments to maintain the same resolution from input file to output file.
# Loop to extract all prism normals and write to tiffs.

extractBox <- function(input_file, lonbds, latbds, input_proj, output_proj, output_res, output_file_path, output_file_name, temp_file_name) {
	require(sp)
	require(rgdal)
	require(raster)
	# make box bounded by lonbds and latbds
	b_box <- SpatialPoints(coords = cbind(lonbds, latbds), proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
	# Project lat long of box into the raster's CRS
	b_box <- spTransform(b_box, CRSobj = CRS(input_proj))

	# define function to make the polygon into a json
	polygon2json <- function(p) {
		paste('{"type": "Polygon", "coordinates": [ [', paste(apply(p, 1, function(x){paste("[",x[1],", ",x[2],"]",sep="")}), collapse=","), '] ]}', sep="")
	}
	
	# Create geojson object from bounding box
	box_json <- polygon2json(as(extent(b_box), 'SpatialPolygons')@polygons[[1]]@Polygons[[1]]@coords)
	box_geojson <- readOGR(box_json, "OGRGeoJSON", verbose = FALSE,  p4s = input_proj)
	# write the geojson to .shp
	writeOGR(box_geojson, output_file_path, temp_file_name, driver="ESRI Shapefile")
	
	# define arguments for system call
	call_args <- paste("-t_srs", output_proj, "-tr", output_res, output_res, "-crop_to_cutline -overwrite -dstnodata NULL -cutline", file.path(output_file_path, paste0(temp_file_name, '.shp')), input_file, file.path(output_file_path,output_file_name))
	# call GDAL to clip the box.
	system2(command="gdalwarp", args=call_args)	
	
	# Remove temporary shapefile
	for (file_ext in c('.shp', '.shx', '.prj', '.dbf')) system2(command="rm", args=file.path(output_file_path, paste0(temp_file_name, file_ext)))
		
}

latbds <- c(41.6, 47.6)
lonbds <- c(-90.5, -82.2)

# Input coordinate system
nad_crs <- '+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0' # for PRISM layers

# Output coordinate system (Michigan GeoRef)
# Note that there are actually quote characters in the below string (unlike the ones above) so that it can be passed as an argument to the system call.
migeoref_crs <- '"+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.255555555556 +k=0.9996 +x_0=2546731.496 +y_0=-4354009.816 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"'

# Loop thru PRISM file names
fp <- '/mnt/research/plz-lab/NEON/external_data/raw_external_data/prismtmp'
vars <- c('ppt','tmax','tmean','tmin')
months <- c(paste0('0',1:9), '10', '11', '12')

vm <- expand.grid(vars, months)

library(foreach)
library(doParallel)
registerDoParallel(cores = 12)

foreach(i = 1:nrow(vm)) %dopar% {
	file_name <- paste('PRISM', vm[i, 1], '30yr_normal_800mM2', vm[i, 2], 'bil', sep = '_')
	extractBox(input_file = file.path(fp, file_name, paste0(file_name, '.bil')),
			   lonbds = lonbds,
			   latbds = latbds,
			   input_proj = nad_crs,
			   output_proj = migeoref_crs,
			   output_res = 30,
			   output_file_path = '/mnt/home/qdr/data/ninarasters',
			   output_file_name = paste('prism', vm[i,1], vm[i,2], 'mi.tif', sep = '_'),
			   temp_file_name = paste0('temp_bbox_', i))
}
