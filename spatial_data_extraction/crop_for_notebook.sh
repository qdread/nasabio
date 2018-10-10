#!/bin/bash
# Create demo rasters for notebook.
# Include a buffer of 100 km around the borders of New Hampshire.

module load GDAL/2.0.1

inpath="/mnt/research/nasabio/data"
outpath="/mnt/research/nasabio/temp/sampledata"
lon1="-73.7"
lon2="-69.5"
lat1="41.5"
lat2="46.5"

# Bioclim 1km resolution (still in lat longs)
gdalwarp -te $lon1 $lat1 $lon2 $lat2 ${inpath}/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.vrt ${outpath}/bioclim1k_NH.tif

# 30m DEM 
gdalwarp -te $lon1 $lat1 $lon2 $lat2 ${inpath}/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big_singlefile.vrt ${outpath}/dem30m_NH.tif

# 100m DEM (in Albers)
gdalwarp -te $lon1 $lat1 $lon2 $lat2 -te_srs '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 100 100 ${inpath}/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt ${outpath}/dem100m_NH.tif

# Geological age (1 km resolution)
gdalwarp -te $lon1 $lat1 $lon2 $lat2 ${inpath}/geology/geo_ages/GEA.vrt ${outpath}/GEA_NH.tif
