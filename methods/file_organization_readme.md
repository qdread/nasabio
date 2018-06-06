# Organization of NASAbioXgeo files on hpcc

QDR, 05 June 2018

This is a description of what files associated with the NASAbioXgeo project are on the MSU server. The files are located at `/mnt/research/nasabio`. I've listed each folder along with a description of what is in it.

### code

This folder has many R scripts and Shell scripts in it. The R scripts are all copies of version-controlled scripts on GitHub that have been uploaded to the server so that the code can be run. It is __*A BAD IDEA*__ to edit R scripts on the HPCC because this can cause versioning conflicts with the GitHub versions. The same goes for the longer Shell scripts in this folder. However some of the shell scripts are `qsub` scripts that are essentially all boilerplate code. The `qsub` scripts are not version-controlled. I have mostly created them directly on the server by copying an older script and modifying a line or two.

### data

Data are kept here. There are a number of subfolders. Most of the subfolders are for the different environmental data layers. The folders `bbs` and `fia` have a lot of files in them, each internally organized with more subdirectories. Most of the files in the subdirectories are temporary files that probably  