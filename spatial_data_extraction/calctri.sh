#!/bin/bash
# Generate TRI, TPI, and roughness from all our geodiversity layers
# Use gdaldem from the command line
# Create TIFs and then use gdalbuildvrt to make a VRT from each one (multi-layer)
# Algorithm option is not supported for TRI, TPI, and roughness.
# Will automatically put NA values in the VRT.

module load GDAL/2.0.1

# Elevation
gdaldem TRI /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt /mnt/research/nasabio/data/tri_tpi/conus_30m_TRI_big.tif -compute_edges
gdaldem TPI /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt /mnt/research/nasabio/data/tri_tpi/conus_30m_TPI_big.tif -compute_edges
gdaldem roughness /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem_big.vrt /mnt/research/nasabio/data/tri_tpi/conus_30m_roughness_big.tif -compute_edges

gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/conus_30m_dem_TRI.vrt /mnt/research/nasabio/data/tri_tpi/conus_30m_TRI_big.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/conus_30m_dem_TPI.vrt /mnt/research/nasabio/data/tri_tpi/conus_30m_TPI_big.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/conus_30m_dem_roughness.vrt /mnt/research/nasabio/data/tri_tpi/conus_30m_roughness_big.tif

# Bioclim 1k (19 bands) 
for i in `seq 1 19`; do
	gdaldem TRI /mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.tif /mnt/research/nasabio/data/tri_tpi/bioclim1k_TRI_${i}.tif -compute_edges -b $i
	gdaldem TPI /mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.tif /mnt/research/nasabio/data/tri_tpi/bioclim1k_TPI_${i}.tif -compute_edges -b $i
	gdaldem roughness /mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.tif /mnt/research/nasabio/data/tri_tpi/bioclim1k_roughness_${i}.tif -compute_edges -b $i
done

gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/bioclim1k_TRI.vrt /mnt/research/nasabio/data/tri_tpi/bioclim1k_TRI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/bioclim1k_TPI.vrt /mnt/research/nasabio/data/tri_tpi/bioclim1k_TPI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/bioclim1k_roughness.vrt /mnt/research/nasabio/data/tri_tpi/bioclim1k_roughness*.tif


# Bioclim5k
for i in `seq 1 19`; do
	gdaldem TRI /mnt/research/nasabio/data/bioclim/Bioclim5k/rasterstack/bioclim5k_20170613.tif /mnt/research/nasabio/data/tri_tpi/bioclim5k_TRI_${i}.tif -compute_edges -b $i
	gdaldem TPI /mnt/research/nasabio/data/bioclim/Bioclim5k/rasterstack/bioclim5k_20170613.tif /mnt/research/nasabio/data/tri_tpi/bioclim5k_TPI_${i}.tif -compute_edges -b $i
	gdaldem roughness /mnt/research/nasabio/data/bioclim/Bioclim5k/rasterstack/bioclim5k_20170613.tif /mnt/research/nasabio/data/tri_tpi/bioclim5k_roughness_${i}.tif -compute_edges -b $i
done

gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/bioclim5k_TRI.vrt /mnt/research/nasabio/data/tri_tpi/bioclim5k_TRI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/bioclim5k_TPI.vrt /mnt/research/nasabio/data/tri_tpi/bioclim5k_TPI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/bioclim5k_roughness.vrt /mnt/research/nasabio/data/tri_tpi/bioclim5k_roughness*.tif

# Biocloud1k
for i in `seq 1 8`; do
	gdaldem TRI /mnt/research/nasabio/data/bioclim/Biocloud1k/biocloud${i}_1k_20170613.tif /mnt/research/nasabio/data/tri_tpi/biocloud1k_TRI_${i}.tif -compute_edges
	gdaldem TPI /mnt/research/nasabio/data/bioclim/Biocloud1k/biocloud${i}_1k_20170613.tif /mnt/research/nasabio/data/tri_tpi/biocloud1k_TPI_${i}.tif -compute_edges
	gdaldem roughness /mnt/research/nasabio/data/bioclim/Biocloud1k/biocloud${i}_1k_20170613.tif /mnt/research/nasabio/data/tri_tpi/biocloud1k_roughness_${i}.tif -compute_edges
done

gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/biocloud1k_TRI.vrt /mnt/research/nasabio/data/tri_tpi/biocloud1k_TRI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/biocloud1k_TPI.vrt /mnt/research/nasabio/data/tri_tpi/biocloud1k_TPI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/biocloud1k_roughness.vrt /mnt/research/nasabio/data/tri_tpi/biocloud1k_roughness*.tif

# Biocloud5k
for i in `seq 1 8`; do
	gdaldem TRI /mnt/research/nasabio/data/bioclim/Biocloud5k/biocloud${i}_5k_20170613.tif /mnt/research/nasabio/data/tri_tpi/biocloud5k_TRI_${i}.tif -compute_edges
	gdaldem TPI /mnt/research/nasabio/data/bioclim/Biocloud5k/biocloud${i}_5k_20170613.tif /mnt/research/nasabio/data/tri_tpi/biocloud5k_TPI_${i}.tif -compute_edges
	gdaldem roughness /mnt/research/nasabio/data/bioclim/Biocloud5k/biocloud${i}_5k_20170613.tif /mnt/research/nasabio/data/tri_tpi/biocloud5k_roughness_${i}.tif -compute_edges
done

gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/biocloud5k_TRI.vrt /mnt/research/nasabio/data/tri_tpi/biocloud5k_TRI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/biocloud5k_TPI.vrt /mnt/research/nasabio/data/tri_tpi/biocloud5k_TPI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/biocloud5k_roughness.vrt /mnt/research/nasabio/data/tri_tpi/biocloud5k_roughness*.tif

# DHI
for i in fpar_v5 gpp_v5 lai8_v5 ndvi_v5; do
	gdaldem TRI /mnt/research/nasabio/data/dhi/dhi_${i}.tif /mnt/research/nasabio/data/tri_tpi/dhi_TRI_${i}.tif -compute_edges
	gdaldem TPI /mnt/research/nasabio/data/dhi/dhi_${i}.tif /mnt/research/nasabio/data/tri_tpi/dhi_TPI_${i}.tif -compute_edges
	gdaldem roughness /mnt/research/nasabio/data/dhi/dhi_${i}.tif /mnt/research/nasabio/data/tri_tpi/dhi_roughness_${i}.tif -compute_edges
done

gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/dhi_TRI.vrt /mnt/research/nasabio/data/tri_tpi/dhi_TRI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/dhi_TPI.vrt /mnt/research/nasabio/data/tri_tpi/dhi_TPI*.tif
gdalbuildvrt -srcnodata -9999 -separate /mnt/research/nasabio/data/tri_tpi/dhi_roughness.vrt /mnt/research/nasabio/data/tri_tpi/dhi_roughness*.tif

# Geology (categorical so no analog I guess, unless we use some categorical neighbor difference sum or something)

# Human Impacts
# footprint
gdaldem TRI /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf.vrt /mnt/research/nasabio/data/tri_tpi/hf_TRI.tif -compute_edges
gdaldem TPI /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf.vrt /mnt/research/nasabio/data/tri_tpi/hf_TPI.tif -compute_edges
gdaldem roughness /mnt/research/nasabio/data/human_impacts/hfp-global-geo-grid/hf.vrt /mnt/research/nasabio/data/tri_tpi/hf_roughness.tif -compute_edges

gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/hf_TRI.vrt /mnt/research/nasabio/data/tri_tpi/hf_TRI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/hf_TPI.vrt /mnt/research/nasabio/data/tri_tpi/hf_TPI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/hf_roughness.vrt /mnt/research/nasabio/data/tri_tpi/hf_roughness.tif

#nightlight
gdaldem TRI /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl.vrt /mnt/research/nasabio/data/tri_tpi/night_TRI.tif -compute_edges
gdaldem TPI /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl.vrt /mnt/research/nasabio/data/tri_tpi/night_TPI.tif -compute_edges
gdaldem roughness /mnt/research/nasabio/data/human_impacts/viirs_nightlights/vcm_orm_ntl.vrt /mnt/research/nasabio/data/tri_tpi/night_roughness.tif -compute_edges

gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/night_TRI.vrt /mnt/research/nasabio/data/tri_tpi/night_TRI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/night_TPI.vrt /mnt/research/nasabio/data/tri_tpi/night_TPI.tif
gdalbuildvrt -srcnodata -9999 /mnt/research/nasabio/data/tri_tpi/night_roughness.vrt /mnt/research/nasabio/data/tri_tpi/night_roughness.tif
