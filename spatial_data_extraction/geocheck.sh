#!/bin/bash
# GEODIVERSITY JOB CHECKER
# Goes through all geo output files, finds the missing ones, and makes a new text file containing qsub commands to run the jobs that are missing (if desired, with a longer walltime)

taxon=$1
walltime=$2

geovars=`awk -F',' '{print $3}' /mnt/research/nasabio/data/geodiv_table_for_gdal.csv`
if [ "$taxon" == "bbs" ]; then
	totaljobs=`awk -F',' '{print $7}' /mnt/research/nasabio/data/geodiv_table_for_gdal.csv`
fi
if [ "$taxon" == "fia" ]; then
	totaljobs=`awk -F',' '{print $8}' /mnt/research/nasabio/data/geodiv_table_for_gdal.csv`
fi

nvars=`wc -l < /mnt/research/nasabio/data/geodiv_table_for_gdal.csv`

cd /mnt/research/nasabio/data/${taxon}/allgeodiv_v2

# Loop through all the variables, then the number of jobs
for j in $(seq 2 $nvars); do
	var=`echo $geovars | awk -v j=$j -F' ' '{print $j}'`
	njob=`echo $totaljobs | awk -v j=$j -F' ' '{print $j}'`
	jobsleft=0
	needtorun=""
	for i in $(seq 1 $njob); do
		# Construct file name and check if it exists.
		# If not, add that job number to the needtorun string.
		fname="${var}_${i}.r"
		if [ ! -e "$fname" ]; then
                ((jobsleft++))
                needtorun+="$i,"
        fi
		# If the number of jobs equals 250 or if we have reached the last job with some to run,
		# create a qsub string and paste it at the bottom of the qsub file.
		if [ "$jobsleft" -eq "250" ] || ([ "$i" -eq "${njob}" ] && [ "$jobsleft" -gt "0" ]); then
			# Create qsub string
			qsub_string="qsub geoextract.sh -N ${taxon}_${var} -v taxon=${taxon},geovar=${var} -l walltime=${walltime},mem=1gb -t ${needtorun%?}"
			# Paste qsub string at the bottom of qsub file
			echo $qsub_string >> /mnt/research/nasabio/code/geo_qsub_all.txt
			# Reset jobs left and needtorun
			jobsleft=0
			needtorun=""
		fi

	done
	
done
