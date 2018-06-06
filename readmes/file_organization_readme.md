# Organization of NASAbioXgeo files on hpcc

QDR, 05 June 2018

This is a description of what files associated with the NASAbioXgeo project are on the MSU server. The files are located at `/mnt/research/nasabio`. I've listed each folder along with a description of what is in it.

### code

This folder has many R scripts and Shell scripts in it. The R scripts are all copies of version-controlled scripts on GitHub that have been uploaded to the server so that the code can be run. It is __*A BAD IDEA*__ to edit R scripts on the HPCC because this can cause versioning conflicts with the GitHub versions. The same goes for the longer Shell scripts in this folder. However some of the shell scripts are `qsub` scripts that are essentially all boilerplate code. The `qsub` scripts are not version-controlled. I have mostly created them directly on the server by copying an older script and modifying a line or two.

There is a subdirectory `routputfiles` where the R console logs are usually written to. I usually clear this folder out every so often. There are a few subdirectories marked `old` with some unused code. Overall, all of the files in the `code` directory are either backed up on GitHub or are small job submission scripts that are easily replaceable. So it's fine to edit or move any file here. 

### data

Data are kept here. There are a number of subfolders. Most of the subfolders are for the different environmental data layers. The folders `bbs` and `fia` have a lot of files in them, each internally organized with more subdirectories. Most of the files in the subdirectories are temporary files that probably won't be needed in the future. However there are a lot of scripts that access them currently, so it's best not to move them around too much or to change the names of the subdirectories.

The organization of the `data` folder is:

#### data/bbs

In this folder there are files related to calculating BBS biodiversity, and geodiversity at BBS points. The raw BBS data are actually stored on the `/mnt/research/aquaxterra` space, so not here. However there are R objects in the directory that have the BBS data saved in a convenient format to access in R (`bbsworkspace_singleyear.r` for by route counts, and `bbsworkspace_bystop_20072016` for by stop counts.) The two folders with clean output CSVs of geodiversity and biodiversity are `geodiversity_CSVs` and `biodiversity_CSVs`. The remaining directories are other raw data sources such as shapefiles, phylogenetic trees, range maps, etc., and temporary files for calculating biodiversity and geodiversity.

#### data/fia

In this folder are all the files related to FIA biodiversity and geodiversity at FIA locations. The raw FIA data (not including true coordinates) are in the folder `treedata10nov`. The true coordinates are not accessible through the group directory. The final geodiversity and biodiversity CSV files are in `geodiversity_CSVs` and `biodiversity_CSVs`. The remaining directories are raw data sources including trait data, phylogenetic trees, and FIA plot condition data by state, or temporary files for calculating biodiversity and geodiversity.

### figs

Figures that need a lot of memory to draw (mostly maps) need to be rendered using the server. They are written here. Most of the figures have been copied to Google Drive.

### restore

This folder was created by the sysadmins back in December 2017 to restore some mistakenly deleted files. It is not needed anymore but it's write-protected so I can't delete it.

### temp

This folder has a lot of temporary files mostly associated with model fitting. There are a number of model fitting scripts I'm working on now that access some of these files, but it is unlikely anyone will need to access them in the future.