#!/bin/sh -login
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -N extrvrtbbs
#PBS -j oe
#PBS -m n
#PBS -t 1-1000

module load R/3.2.0 GDAL GEOS
cd /mnt/research/nasabio/code
now="$(date +'%d%h%Y%H%M')"
R CMD BATCH --no-save --no-restore extractelevfromvrt_bbs.r ./routputfiles/extractelevfromvrt_bbs_${now}.txt
