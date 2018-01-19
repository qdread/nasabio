#!/bin/bash
# BETA-DIVERSITY JOB CHECKER (FIA)
# Goes through all beta output files, finds the missing ones, and appends to the text file containing qsub commands to run the jobs that are missing (if desired, with a longer walltime)

walltime=$1

cd /mnt/research/nasabio/data/fia/diversity/usa

# Loop through all the variables, then the number of jobs
	jobsleft=0
	needtorun=""
	for i in $(seq 1 135174); do
		# Construct file name and check if it exists.
		# If not, add that job number to the needtorun string.
		fname="beta_${i}.r"
		if [ ! -e "$fname" ]; then
                ((jobsleft++))
                needtorun+="$i,"
        fi
		# If the number of jobs equals 250 or if we have reached the last job with some to run,
		# create a qsub string and paste it at the bottom of the qsub file.
		if [ "$jobsleft" -eq "250" ] || ([ "$i" -eq "${njob}" ] && [ "$jobsleft" -gt "0" ]); then
			isbig="yes"
			if [ "$i" -le "100000" ]; then
				isbig="no"
			fi
			# Create qsub string
			qsub_string="qsub fiabd.sh -v isbig=${isbig} -l walltime=${walltime} -t ${needtorun%?}"
			# Paste qsub string at the bottom of qsub file
			echo $qsub_string >> /mnt/research/nasabio/code/bd_qsub.txt
			# Reset jobs left and needtorun
			jobsleft=0
			needtorun=""
		fi

	done
