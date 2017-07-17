#!/bin/bash --login
#PBS -l nodes=1:ppn=1,walltime=00:24:00:00,mem=4gb
#PBS -l feature=gbe
#PBS -j oe
#PBS -N missingmodel
#PBS -m n
#PBS -t 1-249

# Run 3 chains per model, and run the model on each dataset.
DATA_ID=$((PBS_ARRAYID/3 + 1))
CHAIN_ID=$((PBS_ARRAYID%3 + 1))

cd /mnt/home/qdr/code/stancode

module load GNU/6.2

echo "This job is running on $HOSTNAME on `date`"

./phylo_spatial_missing sample num_samples=15000 num_warmup=5000 thin=5 data file=/mnt/research/nasabio/data/fia/staninputs/trait_data_list_${DATA_ID}.R init=0.5 output file=/mnt/research/nasabio/data/fia/stanoutputs/stan_output_${DATA_ID}_${CHAIN_ID}.csv

qstat -f ${PBS_JOBID}
