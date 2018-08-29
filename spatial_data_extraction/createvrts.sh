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
# check out: https://epsg.io/ for specific projections

module load GDAL/2.0.1

cd /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/

###############################################

# MODIS LST: MCMURDO Antarctica: 1km (1000m)
# build multiband vrt out of all singleband annual rasters
gdalbuildvrt -separate /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_MCMURDO.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/*MCMURDO*tif ;  

# reproject and convert vrt to multiband raster
gdalwarp -t_srs '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' -tr 1000 1000 /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_MCMURDO.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_MCMURDO.tif ;

# convert reprojected raster to vrt (necessary?)
gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_MCMURDO.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/mcmurdo/MODIS_MOD11A2_LST_MCMURDO.tif ;

###############################################

# MODIS LST: PALMER Antarctica: 1km (1000m)
# build vrt out of all annual rasters
gdalbuildvrt -separate /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_PALMER.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/*PALMER*tif ;

# reprojecting and converting vrt to multi band tiff
gdalwarp -t_srs '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' -tr 1000 1000 /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_PALMER.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_PALMER.tif ;

# convert reprojected tif to vrt
gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_PALMER.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/palmer/MODIS_MOD11A2_LST_PALMER.tif ;

################################################

# MODIS LST: PUERTO RICO: 1km (1000m)
# clip out Puerto Rico
ls | grep USA | awk -F"_" '{print $5}' | awk -F"." '{print $1}' | while read -r YEAR ; do
	gdalwarp -te -68 -65 17 19 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_USA_"$YEAR".tif /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_PR_"$YEAR".tif ;
done ;

gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_PR.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/*PR*tif ;

gdalwarp -t_srs '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs' -tr 1000 1000 /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_PR.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_PR.tif ;

gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_PR.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_PR.tif ;

###############################################

# MODIS LST: ALASKA: 1km (1000m)
# clip out Alaska
ls | grep USA | awk -F"_" '{print $5}' | awk -F"." '{print $1}' | while read -r YEAR ; do
	gdalwarp -te -157 -139 61 71 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_USA_"$YEAR".tif /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_AK_"$YEAR".tif ;
done ;

gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_AK.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/*AK*tif ;

gdalwarp -t_srs '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs' -tr 1000 1000 /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_AK.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_AK_*.tif ;

gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_AK.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_AK.tif

###############################################

# MODIS LST: USA: 1km (1000m)
gdalbuildvrt -separate /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/*USA*tif ;

gdalwarp -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs' -tr 1000 1000 /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/MODIS_MOD11A2_LST_USA.tif ;

gdalbuildvrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/vrts/MODIS_MOD11A2_LST_USA.vrt /mnt/research/nasabio/data/modis_lst/modis_lst_annual_usa/usa/MODIS_MOD11A2_LST_USA.tif ;
