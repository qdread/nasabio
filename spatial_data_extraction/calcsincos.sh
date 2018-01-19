### run first time only
conda create --yes -n Q gdal python=2.7
python -c "from osgeo import gdal; print(dir(gdal))"
###


source activate Q
cd /mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs

~/code/fia/gdal_calc.py -A conus_30m_aspect_big.vrt --outfile=conus_30m_aspect_sin.tif --calc="sin(A * pi/180)"

~/code/fia/gdal_calc.py -A conus_5k_aspect.vrt --outfile=test.tif --calc="sin(A * pi/180)"

source deactivate Q

###
# on tiles

# sin and cos of each tile
source activate Q
cd /mnt/research/nasabio/data/dem/SRTM_30m_DEM/aspect_tiles

for i in $(seq 1 5); do
	for j in $(seq 1 4); do
		~/code/fia/gdal_calc.py -A aspect_${i}_${j}.tif --outfile=sin_aspect_${i}_${j}.tif --calc="sin(A * pi/180)"
		~/code/fia/gdal_calc.py -A aspect_${i}_${j}.tif --outfile=cos_aspect_${i}_${j}.tif --calc="cos(A * pi/180)"
	done
done

source deactivate Q

# combine tiles into single vrt, then write 1 tif, then write 1 vrt again --- I think this is correct :-O

gdalbuildvrt conus_30m_aspect_sin_manytiles.vrt sin_aspect_*.tif 
gdal_translate conus_30m_aspect_sin_manytiles.vrt conus_30m_aspect_sin.tif
gdalbuildvrt conus_30m_aspect_sin.vrt conus_30m_aspect_sin.tif

gdalbuildvrt conus_30m_aspect_cos_manytiles.vrt cos_aspect_*.tif 
gdal_translate conus_30m_aspect_cos_manytiles.vrt conus_30m_aspect_cos.tif
gdalbuildvrt conus_30m_aspect_cos.vrt conus_30m_aspect_cos.tif
