#!/bin/bash --login
#PBS -l nodes=1:ppn=1,walltime=4:00:00,mem=16gb
#PBS -l feature=gbe
#PBS -j oe
#PBS -N writestansumms
#PBS -m n
#PBS -t 38-160

# parallelize because the summary takes a long time to generate

module load GNU/6.2
cd /mnt/research/nasabio/data/fia/stanoutputs/shorter

readfiles=""

for j in {1..3}
do
read n fname <<< `wc -l stan_output_${PBS_ARRAYID}_${j}.csv`
if [ $n -eq "547" ]; then
	readfiles+="${fname} "
fi
done
if [ ${#readfiles} -gt "0" ]; then
	~/cmdstan-2.15.0/bin/stansummary ${readfiles} > summary_${PBS_ARRAYID}.txt
fi

