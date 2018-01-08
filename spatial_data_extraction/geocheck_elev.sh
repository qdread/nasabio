#!/bin/bash
# GEODIVERSITY JOB CHECKER
# Goes through all geo output files, finds the missing ones, and makes a new text file containing qsub commands to run the jobs that are missing (if desired, with a longer walltime)

walltime=$1

cd /mnt/research/nasabio/data/fia/elevstats_usa

geovars=(elevation slope)
totaljobs=(10000 10000)
mems=(1 1)
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
			qsub_string="qsub elevextract.sh -N ${geovars[j]} -v geovar=${geovars[j]} -l walltime=${walltime},mem=${mems[j]}gb -t ${needtorun%?}"
			# Paste qsub string at the bottom of qsub file
			echo $qsub_string >> /mnt/research/nasabio/code/elev_qsub.txt
			# Reset jobs left and needtorun
			jobsleft=0
			needtorun=""
		fi

	done
	
done
