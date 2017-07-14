#!/bin/bash --login
#PBS -l nodes=1:ppn=1,walltime=00:24:00:00,mem=4gb
#PBS -l feature=gbe
#PBS -j oe
#PBS -N fullmodel_randinit
#PBS -m n
#PBS -t 1-4

cd /mnt/home/f0002182/stan_code

module load GNU/6.2

# Produce a random number between 0.5 and 1.5 for the initial value, using awk since bash only does integers.
initial_value=$(awk "BEGIN {printf \"%.3f\",0.5+(${RANDOM}%1000)/1000}")

echo "This job is running on $HOSTNAME on `date`"

./phylo_spatial_trait sample num_samples=15000 num_warmup=5000 thin=5 data file=./trait_data_list.R init=${initial_value} output file=./stan_output${PBS_ARRAYID}.csv

qstat -f ${PBS_JOBID}
