#!/bin/sh -login
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -N occstan
#PBS -j oe
#PBS -m n
#PBS -t 1-60

CHAIN=$((PBS_ARRAYID%3 + 1))
YEAR=$((PBS_ARRAYID/3 + 1997))

module load GNU/6.2
cd /mnt/research/nasabio/code/occupancy

./occmod_ranef_lesspars sample num_samples=500 num_warmup=1000 thin=1 data file=/mnt/research/nasabio/data/bbs/staninputs/data_ranef${YEAR}.txt init=0.1 output file=/mnt/research/nasabio/data/bbs/stanoutputs/occ_ranef_${YEAR}_${CHAIN}.csv
