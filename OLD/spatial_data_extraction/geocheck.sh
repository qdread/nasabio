#!/bin/bash
# GEODIVERSITY JOB CHECKER
# Goes through all geo output files, finds the missing ones, and makes a new text file containing qsub commands to run the jobs that are missing (if desired, with a longer walltime)

walltime=$1

cd /mnt/research/nasabio/data/fia/allgeodiv

geovars=(bioclim1k bioclim5k biocloud1k biocloud5k dhi elevation slope aspect tpi hf gea night soil)
totaljobs=(30000 6000 30000 6000 1500 135174 135174 135174 135174 1500 1500 1500 1500)
mems=(16 16 16 16 8 80 80 80 80 8 8 8 8)
((nvars=${#geovars[@]}-1))

# Loop through all the variables, then the number of jobs
for j in $(seq 0 $nvars); do
	jobsleft=0
	needtorun=""
	for i in $(seq 1 ${totaljobs[j]}); do
		# Construct file name and check if it exists.
		# If not, add that job number to the needtorun string.
		fname="${geovars[j]}_${i}.r"
		if [ ! -e "$fname" ]; then
                ((jobsleft++))
                needtorun+="$i,"
        fi
		# If the number of jobs equals 250 or if we have reached the last job with some to run,
		# create a qsub string and paste it at the bottom of the qsub file.
		if [ "$jobsleft" -eq "250" ] || ([ "$i" -eq "${totaljobs[j]}" ] && [ "$jobsleft" -gt "0" ]); then
			# Create qsub string
			qsub_string="qsub geoextract.sh -N ${geovars[j]}_fia -v taxon=fia,geovar=${geovars[j]} -l walltime=${walltime},mem=${mems[j]}gb -t ${needtorun%?}"
			# Paste qsub string at the bottom of qsub file
			echo $qsub_string >> /mnt/research/nasabio/code/geo_qsub_all.txt
			# Reset jobs left and needtorun
			jobsleft=0
			needtorun=""
		fi

	done
	
done
