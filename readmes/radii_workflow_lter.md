# Steps for extracting by radius from raster

Created: QDR, 06 June 2018  
Last modified: ACS, 05 Nov 2018

1. Project raster to Albers equal area conic (US/AK/PR) or Polar stereographic (Antarctica), and set desired resolution. Use `gdalwarp` (command-line GDAL 2.0.1 version must be used). If needed refer to `GitHub/nasabio/spatial_data_extraction/LTER/createvrts.sh`
2. Convert to VRT using `gdalbuildvrt` with GDAL from the command line (again must be done with GDAL 2.0.1). If needed refer to `GitHub/nasabio/spatial_data_extraction/LTER/createvrts.sh`
3. Copy VRT, and the TIF it refers to, to the `VRTs_writable` folder on the flash file system: `/mnt/ffs17/groups/nasabio/VRTs_writable`. (note: A script to copy many VRTs in batch is at `GitHub/nasabio/spatial_data_extraction/cpallvrts.sh`)
4. Add rows to lookup table corresponding to your raster, with the file name and resolution set. The version-controlled lookup table is at `GitHub/nasabio/spatial_data_extraction/LTER/geodiv_table_for_gdal_lter.csv`. Copy the edited lookup table to `/mnt/research/nasabio/data`.
5. Create an `sbatch` call to execute `Github/nasabio/spatial_data_extraction/LTER/geoextract_lter.sh`. 
6. Run the `sbatch` call in `/mnt/research/nasabio/code`. This will execute `geoextract_lter.sh` which will in turn execute `master_extract_lter.r` which sources the scripts `extractbox.r` and `extractfromcircle_stack.r`. What is happening is that you are passing arguments to a shell script, which passes arguments to the master R script, which uses GDAL (`gdalwarp --crop-to-cutline`) to create a box around your focal point just larger than the biggest circle you want to extract. GDAL clips the raster to that box and writes a temporary TIF. Then, GDAL creates a series of smaller rasters that are concentric circles around the focal point (or squares with no data outside the circle perimeter). Again a temporary TIF is written for each of the circles. It uses `gdalinfo` to get summary statistics from the rasters. If more complicated statistics beside mean, min, max, and sd are needed, it uses the R package `raster` to load the entire contents of the circles into RAM and calculate fancy stats. 
7. What results is a list of data frames that is saved as a `.r` object. Since the `sbatch` was an array job, each slice is one of the LTER locations (centroids for now). You need to load all the lists into R and concatenate them into a single list. Then bind the rows so that you have a single data frame.
	
