#!/bin/bash
# Code to submit new beta-diversity jobs periodically
# Updated by QDR on 26 Nov 2018 to work with SLURM cluster.

cd /mnt/research/nasabio/code

# Extract number of jobs currently in queue from the output of squeue.
active_jobs=$(squeue -hru $USER | wc -l)

if [ "$active_jobs" -le "500" ]; then
	# Scan the text file until we find a line not containing the word "DONE"
	cmd=`awk '$0 !~ /DONE/ { print; exit }' bd_qsub.txt`
	# Submit that line as a command
	# Edit 26 Jan: test whether the qsub submitted successfully, only edit text file if it did.
	qsubout=$(eval $cmd)
	
	if [ "${#qsubout}" -gt "0" ]; then
		# Add "#DONE" to the end of the line we just submitted
		awk '!f && !/#DONE$/{ $0=$0 " #DONE"; f=1 }1' bd_qsub.txt > tmp && mv -f tmp bd_qsub.txt
	fi
fi

