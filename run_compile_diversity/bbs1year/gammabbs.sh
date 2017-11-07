#!/bin/sh -login
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb
#PBS -N gammabbs_all
#PBS -j oe
#PBS -m n

module load R/3.2.0 GDAL
cp /mnt/home/qdr/code/fia/bbs1year/bbs_allradiigamma.r ${TMPDIR}

cd ${TMPDIR}

R CMD BATCH --no-save --no-restore bbs_allradiigamma.r /mnt/home/qdr/code/fia/routputfiles/bbs_gamma_${PBS_ARRAYID}.txt
