#!/bin/sh -login
#PBS -l walltime=4:00:00	# If you keep this time request 4 hours or less, the job can run on any node so will start faster. Max is 7 days.
#PBS -l mem=4gb				# Keep this low. If it's over about 12 gb, the job will use up more than 1 processor even if it isn't parallel. That will make it start slower.
#PBS -l nodes=1:ppn=1		# This is the number of nodes and processors per node used by each array task. Change ppn only if each single array task will require >1 processor.
#PBS -N jobname				# Change jobname to something informative so you can see it in the queue
#PBS -j oe					# Put error and output into a single file (saves creation of unneccesary files)
#PBS -m n					# Don't send emails
#PBS -t 1-2					# Number of array tasks to be created. Remove this line if it's just a single job.

# Shell script template
# Must include arguments in the -v flag of qsub
# In this template, argument a is a numeric and argument b is a string

module load R/3.2.0 GDAL 	# Here you load all the modules needed for the job. If you need to run R 3.3 instead, comment this line out and uncomment the next 3 lines
# module swap GNU GNU/4.9
# module load OpenMPI/1.10.0
# module load R/3.3.2

cp /mnt/research/nasabio/code/"INSERT NAME OF R SCRIPT HERE" $TMPDIR/rcode.r	# Insert the name of the script here. It gets copied to the temporary job directory
cd $TMPDIR
now="$(date +'%d%h%Y%H%M')"													# Get the current date to be used as the file name for the R console log

# Paste together the arguments and execute the pasted string as a command.
# This will write a console log for each array job to the routputfiles directory. Often helpful to refer to it to see if errors came up.
cmd="R CMD BATCH --no-save --no-restore '--args a="$a" b=\""$b"\"' rcode.r /mnt/research/nasabio/code/routputfiles/jobname_"$PBS_ARRAYID"_"$now".txt"
eval $cmd

##### QSUB FOR THIS JOB
# In the terminal, cd to the directory where this script is located and run:
# qsub job_script_template.sh -v a=50,b=yes
# That will execute the R script with arguments a=50 and b='yes' as long as you parse the arguments in the head of the R script
# You can also change other arguments in the qsub call to override the ones at the top of the script.
# For example if you wanted to run it with different arguments, higher memory usage, and only run task number 7 you could call
# qsub job_script_template.sh -v a=500,b=no -l mem=8gb -t 7

# Also see https://wiki.hpcc.msu.edu/display/hpccdocs/HPCC+Quick+Reference+Sheet

