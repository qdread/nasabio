module load GDAL/2.0.1

gdalwarp -t_srs '+proj=utm +zone=18 +south +ellps=intl +towgs84=-288,175,-376,0,0,0,0 +units=m +no_defs' -tr 1000 1000 -te -78 -5 -71 9.5 -te_srs '+proj=longlat +ellps=WGS84 +no_defs' /mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.vrt /mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_colombia_ecuador.tif