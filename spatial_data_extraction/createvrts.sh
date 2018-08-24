#!/bin/bash
# TITLE:         Prepare environmental VRTS for LTER sites (vrts enable working with very large rasters)
# AUTHOR:        QDR / modified by PLZ/ ACS 
# COLLABORATORS: Sydne Record, Eric Sokol, Riley Andrade, ... LTER metacommunities synth grp
# DATA:          LTER site to regional level environmental data 
# PROJECT:       "LTER Metacommunities Synthesis Working Group"
# DATE:          initiated: August 2018; last run:

# Use Albers Equal Area Conic.
# GDAL from command line

# Raster data: mean MODIS LST MOD11A2; max NDVI MOD13Q1; 
# downloaded from Google Earth Engine 21-22 August 2018
# Sites span AK to Antarctica so they need different projections
# check out: https://epsg.io/3031 for specific projections

module load GDAL/2.0.1

# MODIS LST: MCMURDO Antarctica: 1km (1000m)
gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_MCMURDO.vrt *.tif 

gdalwarp -t_srs '+proj=aea +lat_1=-71 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 1000 1000 /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_MCMURDO.vrt MODIS_MOD11A2_LST_MCMURDO.tif

gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_MCMURDO.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_MCMURDO.tif

# MODIS LST: USA: 1km (1000m)
gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA_${i}.tif

gdalwarp -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 1000 1000 /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA_*.tif

gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA.tif



# ones with global extent don't work for Albers reprojection.
# must first crop to cutline of entire USA.
# For our "big" USA map, it's  -te -125 20 -65 56 (xmin ymin xmax ymax)
# Must specify that target extent is given in lat longs, though the output file will be in albers. (use -te_srs flag)

# DHI (four bands, each to a tif, then to a single vrt)
for i in fpar gpp lai8 ndvi; do
	gdalwarp -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 5000 5000 -te -125 20 -65 56 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/dhi/dhi_${i}.vrt /mnt/research/nasabio/data/dhi/dhi_5k_${i}.tif 
done

gdalbuildvrt -separate -b 1 /mnt/research/nasabio/data/dhi/dhi_5k.vrt /mnt/research/nasabio/data/dhi/dhi_5k_*.tif 

# Human footprint
gdalwarp -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 5000 5000 -te -125 20 -65 56 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf.vrt /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf_5k.tif

gdalbuildvrt /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf_5k.vrt /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf_5k.tif

# Nightlight
gdalwarp -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 5000 5000 -te -125 20 -65 56 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl.vrt /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl_5k.tif

gdalbuildvrt /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl_5k.vrt /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl_5k.tif

#####################################################
# Generate TRI, TPI, and roughness metrics for the coarser rasters.

# Elevation
gdaldem TRI /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_5k_dem.vrt /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_TRI.tif -compute_edges
gdaldem TPI /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_5k_dem.vrt /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_TPI.tif -compute_edges
gdaldem roughness /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_5k_dem.vrt /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_roughness.tif -compute_edges

gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_TRI.vrt /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_TRI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_TPI.vrt /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_TPI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_roughness.vrt /mnt/research/nasabio/data/tri_tpi/conus_5k_dem_roughness.tif

# DHI
for i in fpar gpp lai8 ndvi; do
	gdaldem TRI /mnt/research/nasabio/data/dhi/dhi_5k_${i}.tif /mnt/research/nasabio/data/tri_tpi/dhi_5k_TRI_${i}.tif -compute_edges
	gdaldem TPI /mnt/research/nasabio/data/dhi/dhi_5k_${i}.tif /mnt/research/nasabio/data/tri_tpi/dhi_5k_TPI_${i}.tif -compute_edges
	gdaldem roughness /mnt/research/nasabio/data/dhi/dhi_5k_${i}.tif /mnt/research/nasabio/data/tri_tpi/dhi_5k_roughness_${i}.tif -compute_edges
done

gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/dhi_5k_TRI.vrt /mnt/research/nasabio/data/tri_tpi/dhi_5k_TRI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/dhi_5k_TPI.vrt /mnt/research/nasabio/data/tri_tpi/dhi_5k_TPI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/dhi_5k_roughness.vrt /mnt/research/nasabio/data/tri_tpi/dhi_5k_roughness*.tif

# Human footprint
gdaldem TRI /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf_5k.vrt /mnt/research/nasabio/data/tri_tpi/hf_5k_TRI.tif -compute_edges
gdaldem TPI /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf_5k.vrt /mnt/research/nasabio/data/tri_tpi/hf_5k_TPI.tif -compute_edges
gdaldem roughness /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf_5k.vrt /mnt/research/nasabio/data/tri_tpi/hf_5k_roughness.tif -compute_edges

gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/hf_5k_TRI.vrt /mnt/research/nasabio/data/tri_tpi/hf_5k_TRI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/hf_5k_TPI.vrt /mnt/research/nasabio/data/tri_tpi/hf_5k_TPI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/hf_5k_roughness.vrt /mnt/research/nasabio/data/tri_tpi/hf_5k_roughness.tif

# Nightlight
gdaldem TRI /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl_5k.vrt /mnt/research/nasabio/data/tri_tpi/night_5k_TRI.tif -compute_edges
gdaldem TPI /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl_5k.vrt /mnt/research/nasabio/data/tri_tpi/night_5k_TPI.tif -compute_edges
gdaldem roughness /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl_5k.vrt /mnt/research/nasabio/data/tri_tpi/night_5k_roughness.tif -compute_edges

gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/night_5k_TRI.vrt /mnt/research/nasabio/data/tri_tpi/night_5k_TRI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/night_5k_TPI.vrt /mnt/research/nasabio/data/tri_tpi/night_5k_TPI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/night_5k_roughness.vrt /mnt/research/nasabio/data/tri_tpi/night_5k_roughness.tif

######################################################
# Added 10 Jan: use mode to resample the Geological Age and Soil Type rasters to 5 km.

gdalwarp -r mode -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 5000 5000 -te -125 20 -65 56 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/geology/geo_ages/GEA.vrt /mnt/research/nasabio/data/geology/geo_ages/GEA_5k.tif

gdalbuildvrt /mnt/research/nasabio/data/geology/geo_ages/GEA_5k.vrt /mnt/research/nasabio/data/geology/geo_ages/GEA_5k.tif

gdalwarp -r mode -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0' -tr 5000 5000 -te -125 20 -65 56 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/geology/soils/stg.vrt /mnt/research/nasabio/data/geology/soils/stg_5k.tif

gdalbuildvrt /mnt/research/nasabio/data/geology/soils/stg_5k.vrt /mnt/research/nasabio/data/geology/soils/stg_5k.tif

#########################################################
# Added 06 June: convert Idaho fire to VRT
# No need to reproject b/c already in Albers
gdalbuildvrt /mnt/research/nasabio/data/mtbs_ID_2012.vrt /mnt/research/nasabio/data/mtbs_ID_2012.tif