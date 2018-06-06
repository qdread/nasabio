# HPCC tricks of the trade

QDR, 06 June 2018

## File systems and locations

These file systems relevant to our our work, with "do's and don'ts" for each one.

### Research space

The group's research space is `/mnt/research/nasabio`. It has 2 TB capacity, as of June 2018 about 1 TB is full.

**DO** store raw data, scripts, temporary data-processing files (.RData and .csv) and final output data (.csv)  
**DO** run R scripts that read a file from this space at the beginning, do something, then save the output. For example, a model fitting script that loads data, fits the model, and saves the results.   
**DON'T** edit any scripts or lookup tables here. Any editing should be done in the GitHub repository where we can maintain version control, and then the file should be uploaded to the HPCC with SFTP. That will keep any conflicts to a minimum.  
**DON'T** run many R scripts in parallel that read one or more files in this space many times and create lots of temporary files in this space. This destabilizes the entire system and brings the wrath of sysadmins on us.

### Scratch space

The group's scratch space is `/mnt/ls15/scratch/groups/nasabio`. It is only accessible from a dev node.

**DO** Store huge temporary files here since the storage space is effectively unlimited.  
**DON'T** Keep anything here that you need in the future that isn't backed up, since unused files disappear after 45 days.  
**DON'T** run many R scripts in parallel that read one or more files in this space many times and create lots of temporary files in this space. This destabilizes the entire system and brings the wrath of sysadmins on us.

### Flash file system

The group's flash file system is `/mnt/ffs17/groups/nasabio`. It has 500 GB capacity.

**DO** copy any files you need to read in parallel lots of times to this directory. You can run R scripts that read the files in parallel over and over without destabilizing the system.    
**DON'T** Keep anything permanent here that isn't backed up. Because of the relatively smaller size of the directory we have to delete unneeded files from it to make room for new ones.

### Temporary job space

There is a special directory &tilde;400 GB in size that is created whenever a job is run. If the job is an array job, a different directory is made for each task. The directory has a different name every time but its name is stored in the system variable `$TMPDIR`.

**DO** Write temporary files, including rasters, to this folder in the course of your job, as long as they are not being read again in parallel lots of times.  
**DO** Copy the R (or other language) script being run by your job script to this folder, because in case any large temporary R environment file is created, it will just disappear when the job is done if the R script is called from this folder. This is a good idea to cover your butt--it will keep a lot of huge hidden files from being created in the research space and unexpectedly filling up the quota.    
**DON'T** Get too attached to anything in here, since it disappears when the job is done, like a daemon from *The Golden Compass*. You have to save any output from the job to the research space or scratch space.

## File permissions

For some reason, the default permissions when I create a new file are for me to have read and write access, the group to have read-only access, and all users to have read-only access. (`rw-r--r--`, or 644 using the permission codes). I would recommend setting it so anyone in our group can read and write, and no one else (`rw-rw----` or 660). That would be done using `chmod 660 <filename>` or `chmod -R 660 <directoryname>`. If you want the file to be executable too, use 770 instead of 660. 

## Issues with job submission

These are the different tricks that I figured out to get jobs to run as quickly and efficiently as possible.

- *Split jobs into tasks that can run in <4 hours.* We have priority access to 20 processors, which we paid for. The maximum job length is 7 days so we can get 20 jobs started right away that are up to 7 days long. However we often need to use more than 20. Communal nodes are available based on priority, but you can only get on other people's private nodes if your job is 4 hours long or less. That's why it's a good idea to run array jobs that have a number of tasks, each of which can run in 4 hours or less. The tasks will start and finish more quickly than one long job. 
- *Use array jobs for parallelizing when possible, not R packages like doParallel*. While the R parallel packages are great, it requires you to submit a job that requests a lot of processors. It's quicker if you just run the same R script with many copies using an array job. That's because you are basically submitting a lot of jobs each requesting only 1 processor. I still use `doParallel` for some purposes but use it sparingly.
- *Minimize memory demands*. Your jobs will start more quickly if you request less memory. Anything above 40 GB or so per task starts to make the tasks slower to start, because you are starting to require multiple processors to get that much RAM, even if the job isn't parallel. Sometimes it's unavoidable but you should keep it in mind. The best way to do this is to not use R because it is such a hog of RAM.
- *Make sure the correct directory is used every time*. The job submission script I made to run R scripts starts by copying the R script into the temporary directory so that the R environment will be there, specifies the directory where the R console log will be written to. Any files that have to be read a lot of times have to be in the flash file system, any temporary files are written to the temporary directory, and the final result is written to the research space.

## Passing arguments from qsub to shell to R

This is kind of tricky but you basically have to include arguments, unquoted, in the `-v` flag of the `qsub` call. Then use the `--args` flag of the `R CMD BATCH` call to pass those arguments to the R script. You need to use string concatenation in Bash to pass those arguments, because the `--args` section has to have single quotes around it, and the individual arguments have to have quotation marks around them if they're strings. Then the R script needs to have some boilerplate code at the top to parse those arguments. See `https://www.r-bloggers.com/including-arguments-in-r-cmd-batch-mode/`.

## Automating job submission with crontab

It's not possible to just dump all the `qsub` calls to the terminal at once because of 2 limitations: each array job can only have 250 tasks, and a user can only have 1000 tasks in their queue. That limits you to only submit 4 max-sized array jobs at once. If you don't have too many jobs, you can manually submit them whenever your job queue clears out. However, if there are a lot, I made some shell scripts to automate the process.

- Put all the `qsub` calls into a text file, one per line, and put the file in `mnt/research/nasabio/code`.
- Edit the `autogeo.sh` script in the `code` folder to refer to the correct text file.
- Set a `crontab` on one of the dev nodes to execute `autogeo.sh` at regular intervals, for example every 30 minutes. Whenever `autogeo.sh` runs, it checks your job queue to see if less than 750 jobs are in the queue. If so, it reads the text file of `qsub` calls until it finds one that does not end with `DONE`. It submits that line and then appends `DONE` to the line. 