# TIFs to single VRT
# GDAL from command line

module load GDAL/2.0.1

cd /mnt/research/nasabio/data/bioclim/Biocloud5k
gdalbuildvrt -srcnodata -Inf -vrtnodata NULL -separate biocloud5k.vrt *.tif

cd /mnt/research/nasabio/data/bioclim/Biocloud1k
gdalbuildvrt -srcnodata -Inf -vrtnodata NULL -separate biocloud1k.vrt *.tif

cd /mnt/research/nasabio/data/dhi
gdalbuildvrt -srcnodata -Inf -vrtnodata NULL -separate dhi.vrt *.tif