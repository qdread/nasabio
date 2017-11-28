# Test extract box transforming to AEA with 30 meter resolution.

lon1 <- -98.59595
lat1 <- 39.79797

lon2 <- -120
lat2 <- 37.5

elevfile <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt'

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fp <- '/mnt/research/nasabio/temp'

extractBox(lon=lon1,
		   lat=lat1,
		   radius=300,
		   input_file=elevfile,
		   input_proj=wgs_crs,
		   output_proj=aea_crs,
		   output_res=30,
		   output_file_path=fp,
		   output_file_name='testbboxaea1.tif')

extractBox(lon=lon2,
		   lat=lat2,
		   radius=300,
		   input_file=elevfile,
		   input_proj=wgs_crs,
		   output_proj=aea_crs,
		   output_res=30,
		   output_file_path=fp,
		   output_file_name='testbboxaea2.tif')		   
		   
x1 <- raster(file.path(fp, 'testbboxaea1.tif'))
x2 <- raster(file.path(fp, 'testbboxaea2.tif'))