# Extract pixel values from different-sized circlees around focal coordinates.
# Create temporary copy of VRT to try to parallelize without crashing.

n_slices <- 1000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# FIA lat long coordinates
library(dplyr)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

# Function to do the extracting
source('/mnt/research/nasabio/code/extractfromcircle.r')

copy_vrt <- function(old_file, new_file) system2("cp", args = paste(old_file, new_file))

vrt_dir <- '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs'
old_vrt <- 'conus_30m_dem.vrt'
new_vrt <- paste0('fia_elev_vrt_', slice, '.vrt')
copy_vrt(old_file = file.path(vrt_dir, old_vrt), new_file = file.path(vrt_dir, new_vrt))

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Radii of circles where we want summary stats
radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300) # in km

# Call extraction function, specifying the raster from which to extract data.
stats_by_point <- extractFromCircle(coords = with(fiacoords, cbind(lon, lat))[rowidxmin:rowidxmax,],
									raster_file = file.path(vrt_dir, new_vrt),
									radii = radii,
									fp = '/mnt/research/nasabio/temp',
									filetag = paste0('fia',slice))

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats/30m/stats_',slice,'.r'))
system2('rm', args = file.path(vrt_dir, new_vrt))
