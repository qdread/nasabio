# GDAL operations 1/8/2018

module load GDAL/2.0.1
cd /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs

# Create VRTs that are small relative to the TIFs
gdal_translate conus_30m_slope_big.vrt conus_30m_slope_big.tif
gdal_translate conus_30m_aspect_big.vrt conus_30m_aspect_big.tif
gdalbuildvrt -overwrite conus_30m_slope_big.vrt conus_30m_slope_big.tif
gdalbuildvrt -overwrite conus_30m_aspect_big.vrt conus_30m_aspect_big.tif

gdal_translate conus_30m_dem_big.vrt conus_30m_dem_big.tif
gdalbuildvrt conus_30m_dem_big_singlefile.vrt conus_30m_dem_big.tif

# Calculate slope and aspect at 5km resolution
gdaldem slope conus_5k_dem.vrt conus_5k_slope.tif
gdaldem aspect conus_5k_dem.vrt conus_5k_aspect.tif
gdalbuildvrt conus_5k_slope.vrt conus_5k_slope.tif
gdalbuildvrt conus_5k_aspect.vrt conus_5k_aspect.tif