# HDF test

# BBS lat long coordinates
bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')

raster_file = '/mnt/research/nasabio/data/bioclim/Bioclim5k/rasterstack/bioclim5k_20170613.tif'
raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt'
radii = 50
r = 1
fp = '/mnt/research/nasabio/temp'
filetag = 'bbstest'
nlayers = 19

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

circle_r <- latLongCircle(loni, lati, radii[r], 1000)
circle_json <- polygon2json(circle_r)
circle_geojson <- readOGR(circle_json, "OGRGeoJSON", verbose = F,  p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# write the geojson to .shp
writeOGR(circle_geojson, fp, paste0("temp_circle_", filetag), driver="ESRI Shapefile")
# define arguments for system call
call_args <- paste("-crop_to_cutline -overwrite -dstnodata NULL -cutline", file.path(fp, paste0("temp_circle_", filetag, ".shp")), raster_file, file.path(fp,paste0("temp_circle_extracted", filetag, ".hdf")))
system2(command="gdalwarp", args=call_args)
g.info <- system2(command="gdalinfo", args=paste("-stats", file.path(fp,paste0("temp_circle_extracted", filetag, ".hdf"))), stdout=TRUE) 
g.info
system2(command="rm", args=file.path(fp, paste0("temp_circle*", filetag, "*")))