#!/bin/bash
# FULLY AUTOMATED GEODIVERSITY CALCULATION
# Runs every 20 minutes for a maximum of 72*250 jobs that can be submitted daily by one account.
# If total number of jobs <= 750, it reads the .txt file containing the list of qsubs until it finds the first line not containing #DONE
# It submits that line and then appends #DONE to the line in that .txt file so that another account can't resubmit the same one.

cd /mnt/research/nasabio/code

# Extract number of jobs currently in queue from the output of showq.
active_jobs=`/opt/moab/bin/showq -u $USER | awk '/^Total/ { print $3 }'`

if [ "$active_jobs" -le "750" ]; then
	# Scan the text file until we find a line not containing the word "DONE"
	cmd=`awk '$0 !~ /DONE/ { print; exit }' geo_qsub_all.txt`
	# Submit that line as a command
	# Edit 26 Jan: test whether the qsub submitted successfully, only edit text file if it did.
	qsubout=$(eval $cmd)
	# Add "#DONE" to the end of the line we just submitted
	if [ "${#qsubout}" -gt "0" ]; then
		awk '!f && !/#DONE$/{ $0=$0 " #DONE"; f=1 }1' geo_qsub_all.txt > tmp && mv -f tmp geo_qsub_all.txt
		chmod a+w geo_qsub_all.txt
	fi
fi
