#!/bin/bash --login
#PBS -l nodes=1:ppn=1,walltime=4:00:00,mem=16gb
#PBS -l feature=gbe
#PBS -j oe
#PBS -N writestansumms
#PBS -m n
#PBS -t 1-83

# parallelize because the summary takes a long time to generate

module load GNU/6.2
cd /mnt/research/nasabio/data/fia/stanoutputs/shorter

# use wc -l to see if all three chains have 547 lines

read n1 f1 n2 f2 n3 f3 <<< `wc -l stan_output_${PBS_ARRAYID}_*.csv`

if [ $n1 -eq "547" -a $n2 -eq "547" -a $n3 -eq "547" ]; then
	~/cmdstan-2.15.0/bin/stansummary stan_output_${PBS_ARRAYID}_*.csv > summary_${PBS_ARRAYID}.txt
fi