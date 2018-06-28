### COPY ALL TIF AND VRT TO FLASH FILE SYSTEM DIRECTORY
### (Only 500GB can be stored there at any given time)

TO="/mnt/ffs17/groups/nasabio/VRTs"
FROM="/mnt/research/nasabio/data"

# Bioclim and Biocloud
cp ${FROM}/bioclim/Bioclim1k/rasterstack/* ${TO}
cp ${FROM}/bioclim/Bioclim5k/rasterstack/* ${TO}
cp ${FROM}/bioclim/Biocloud1k/* ${TO}
cp ${FROM}/bioclim/Biocloud5k/* ${TO}

# DEM
cp ${FROM}/dem/SRTM_30m_DEM/VRTs/* ${TO}
cp ${FROM}/dem/SRTM_30m_DEM/aspect_tiles/*sin* ${TO}
cp ${FROM}/dem/SRTM_30m_DEM/aspect_tiles/*cos* ${TO}

# DHI
cp ${FROM}/dhi/dhi* ${TO}

# Geology
cp ${FROM}/geology/geo_ages/GEA_5k* ${TO}
cp ${FROM}/geology/geo_ages/GEAISG3a.tif ${TO}
cp ${FROM}/geology/geo_ages/GEA.vrt ${TO}
cp ${FROM}/geology/soils/stg_5k* ${TO}
cp ${FROM}/geology/soils/STGHWS1a.tif ${TO}
cp ${FROM}/geology/soils/stg.vrt ${TO}

# Humans
cp ${FROM}/human_impacts/hfp-global-geo-grid/hf_5k* ${TO}
cp ${FROM}/human_impacts/hfp-global-geo-grid/hf_v2geo1.tif ${TO}
cp ${FROM}/human_impacts/hfp-global-geo-grid/hf.vrt ${TO}
cp ${FROM}/human_impacts/viirs_nightlights/vcm* ${TO}
cp ${FROM}/human_impacts/viirs_nightlights/SVDNB_npp_20150101-20151231_75N180W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif ${TO}

# TRI (lots)
cp ${FROM}/tri_tpi/*roughness* ${TO}
cp ${FROM}/tri_tpi/*TRI* ${TO}
