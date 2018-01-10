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

# Use anaconda to do raster arithmetic 10 Jan 2018

source ~/.bashrc

try:
    from osgeo import gdal
except ImportError:
    import gdal
try:
    from osgeo import gdalnumeric
except ImportError:
    import gdalnumeric

export PATH="/mnt/home/qdr/anaconda2/bin:/opt/software/Stata/15.0/SE:/opt/software/R/2.15.1--GCC-4.4.5/bin:/opt/software/MATLAB/R2014a/bin:/opt/software/cmake/2.8.5--GCC-4.4.5/bin:/opt/software/OpenMPI/1.4.3--GCC-4.4.5/bin:/usr/lib64/qt-3.3/bin:/opt/software/lmod/bin:/opt/moab/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/local/hpcc/bin:/opt/ibutils/bin"
	
	
cd /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs

~/code/fia/gdal_calc.py -A conus_5k_aspect.tif --outfile=conus_5k_aspect_cos.tif --calc="cos(A*pi/180)"
