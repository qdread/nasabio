#!/bin/bash
# GEODIVERSITY JOB CHECKER
# Goes through all geo output files, finds the missing ones, and makes a new text file containing qsub commands to run the jobs that are missing (if desired, with a longer walltime)

taxon=$1
walltime=$2

cd /mnt/research/nasabio/data/${taxon}/allgeodiv_v2

geovars=(bioclim1k bioclim5k biocloud1k biocloud5k dhi elevation slope hf night bioclim1k_tri bioclim1k_tpi bioclim1k_roughness bioclim5k_tri bioclim5k_roughness biocloud1k_tri biocloud1k_roughness biocloud5k_tri biocloud5k_roughness tri roughness tri_5k roughness_5k aspect_sin_5k aspect_cos_5k slope_5k dhi_tri dhi_roughness dhi_5k_tri dhi_5k_roughness hf_tri hf_roughness hf_5k_tri hf_5k_roughness night_tri night_roughness night_5k_tri night_5k_roughness elevation_5k dhi_5k hf_5k night_5k aspect_sin aspect_cos gea soil gea_5k)
if [ "$taxon" -e "bbs" ]; then
	totaljobs=(500 250 500 250 100 500 500 100 100 500 500 500 250 250 500 500 250 250 500 500 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 500 500 100 100 100)
fi
if [ "$taxon" -e "fia" ]; then
	totaljobs=(5000 1000 5000 1000 500 10000 10000 500 500 5000 5000 5000 1000 1000 5000 5000 1000 1000 10000 10000 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 10000 10000 500 500 500)
fi
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
			qsub_string="qsub geoextract.sh -N ${taxon}_${geovars[j]} -v taxon=${taxon},geovar=${geovars[j]} -l walltime=${walltime},mem=1gb -t ${needtorun%?}"
			# Paste qsub string at the bottom of qsub file
			echo $qsub_string >> /mnt/research/nasabio/code/geo_qsub_all.txt
			# Reset jobs left and needtorun
			jobsleft=0
			needtorun=""
		fi

	done
	
done
