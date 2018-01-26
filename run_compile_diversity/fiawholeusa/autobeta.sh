#!/bin/bash
# Code to submit new beta-diversity jobs every 60 minutes

cd /mnt/research/nasabio/code

# Extract number of jobs currently in queue from the output of showq.
active_jobs=`/opt/moab/bin/showq -u $USER | awk '/^Total/ { print $3 }'`

if [ "$active_jobs" -le "750" ]; then
        # Scan the text file until we find a line not containing the word "DONE"
        cmd=`awk '$0 !~ /DONE/ { print; exit }' bd_qsub.txt`
        # Submit that line as a command
        # Edit 26 Jan: test whether the qsub submitted successfully, only edit text file if it did.
		qsubout=$(eval $cmd)
        # Add "#DONE" to the end of the line we just submitted
		if [ "${#qsubout}" -gt "0" ]; then
			awk '!f && !/#DONE$/{ $0=$0 " #DONE"; f=1 }1' bd_qsub.txt > tmp && mv -f tmp bd_qsub.txt
			chmod a+w bd_qsub.txt
		fi
fi

