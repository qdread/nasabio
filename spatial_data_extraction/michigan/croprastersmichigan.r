# Crop Michigan from different rasters and save the results as new, smaller rasters.
# Edit: added arguments to maintain the same resolution from input file to output file.

extractBox <- function(input_file, lonbds, latbds, input_proj, output_proj, output_res, output_file_path, output_file_name) {
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
	writeOGR(box_geojson, output_file_path, "temp_bbox", driver="ESRI Shapefile")
	
	# define arguments for system call
	call_args <- paste("-t_srs", output_proj, "-tr", output_res, output_res, "-crop_to_cutline -overwrite -dstnodata NULL -cutline", file.path(output_file_path, "temp_bbox.shp"), input_file, file.path(output_file_path,output_file_name))
	# call GDAL to clip the box.
	system2(command="gdalwarp", args=call_args)	
	
	# Remove temporary shapefile
	for (file_ext in c('.shp', '.shx', '.prj', '.dbf')) system2(command="rm", args=file.path(output_file_path, paste0("temp_bbox", file_ext)))
		
}

latbds <- c(41.6, 47.6)
lonbds <- c(-90.5, -82.2)

# Input coordinate systems (varies by dataset)
wgs_crs <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' # for elevation layers
nad_crs <- '+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0' # for PRISM layers
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0' # for NLCD

# Output coordinate system (Michigan GeoRef)
# Note that there are actually quote characters in the below string (unlike the ones above) so that it can be passed as an argument to the system call.
migeoref_crs <- '"+proj=omerc +lat_0=45.30916666666666 +lonc=-86 +alpha=337.255555555556 +k=0.9996 +x_0=2546731.496 +y_0=-4354009.816 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"'

# Extract NLCD 2011
extractBox(input_file = '/mnt/research/aquaxterra/DATA/reprojected_data/NLCD/nlcd_2011.tif',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = aea_crs,
		   output_proj = migeoref_crs,
		   output_res = 30,
		   output_file_path = '/mnt/home/qdr/data/ninarasters',
		   output_file_name = 'nlcd2011_mi.tif')

# Extract aspect DEM		   
extractBox(input_file = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_aspect_big.vrt',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = wgs_crs,
		   output_proj = migeoref_crs,
		   output_res = 30,
		   output_file_path = '/mnt/home/qdr/data/ninarasters',
		   output_file_name = 'aspect_30m_mi.tif')

# Extract slope DEM		   
extractBox(input_file = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_slope_big.vrt',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = wgs_crs,
		   output_proj = migeoref_crs,
		   output_res = 30,
		   output_file_path = '/mnt/home/qdr/data/ninarasters',
		   output_file_name = 'slope_30m_mi.tif')
		   
# Extract elevation DEM		   
extractBox(input_file = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = wgs_crs,
		   output_proj = migeoref_crs,
		   output_res = 30,
		   output_file_path = '/mnt/research/plz-lab/DATA/HWA_MISGP/statewide_hemlock_map/raw_data',
		   output_file_name = 'elevation_30m_mi.tif')

extractBox_noreproject <- function(input_file, lonbds, latbds, input_proj, output_file_path, output_file_name, temp_file_name = 'temp_bbox') {
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
	call_args <- paste("-crop_to_cutline -overwrite -dstnodata NULL -cutline", file.path(output_file_path, paste0(temp_file_name, '.shp')), input_file, file.path(output_file_path,output_file_name))
	# call GDAL to clip the box.
	system2(command="gdalwarp", args=call_args)	
	
	# Remove temporary shapefile
	for (file_ext in c('.shp', '.shx', '.prj', '.dbf')) system2(command="rm", args=file.path(output_file_path, paste0(temp_file_name, file_ext)))
		
}
		   
extractBox_noreproject(input_file = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = wgs_crs,
		   output_file_path = '/mnt/research/plz-lab/DATA/HWA_MISGP/statewide_hemlock_map/raw_data',
		   output_file_name = 'elevation_30m_mi_latlong.tif')

extractBox_noreproject(input_file = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_aspect_big.vrt',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = wgs_crs,
		   output_file_path = '/mnt/research/plz-lab/DATA/HWA_MISGP/statewide_hemlock_map/raw_data',
		   output_file_name = 'aspect_30m_mi_latlong.tif')

extractBox_noreproject(input_file = '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_slope_big.vrt',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = wgs_crs,
		   output_file_path = '/mnt/research/plz-lab/DATA/HWA_MISGP/statewide_hemlock_map/raw_data',
		   output_file_name = 'slope_30m_mi_latlong.tif')
		   
# Test one prism layer
extractBox(input_file = '/mnt/research/plz-lab/NEON/external_data/raw_external_data/prismtmp/PRISM_ppt_30yr_normal_800mM2_01_bil/PRISM_ppt_30yr_normal_800mM2_01_bil.bil',
		   lonbds = lonbds,
		   latbds = latbds,
		   input_proj = nad_crs,
		   output_proj = migeoref_crs,
		   output_res = 30,
		   output_file_path = '/mnt/home/qdr/data/ninarasters',
		   output_file_name = 'prism_ppt_01_mi.tif')


