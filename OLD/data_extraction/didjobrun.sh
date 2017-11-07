#!/bin/bash

direc=$1
prefix=$2
suffix=$3
totaljobs=$4
jobscriptfile=$5

cd $direc

jobsleft=0
needtorun=""

# Loop through and find all the jobs that didn't run.
for i in $(seq 1 $totaljobs); do
	fname="$prefix$i$suffix"
	if [ ! -e "$fname" ]; then
		((jobsleft++))
		needtorun+="$i,"
	fi
done	

echo "$jobsleft jobs out of a total of $totaljobs still need to run."
echo "${needtorun%?}"

# Edit the job submission script, replacing the line with job numbers.
cat $jobscriptfile | sed "s/PBS -t.*/PBS -t ${needtorun%?}/g" > ./dummy.sh
mv dummy.sh $jobscriptfile

# Submit the new jobs!
qsub $jobscriptfile
