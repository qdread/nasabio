# Extract pixel values from different-sized circlees around focal coordinates.
# 250 parallel jobs to extract 30-m srtm data for fia plot coordinates.

n_slices <- 250 
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

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

# Radii of circles where we want summary stats
radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300, 400, 500) # in km

# Call extraction function, specifying the raster from which to extract data.
stats_by_point <- extractFromCircle(coords = with(fiacoords, cbind(lon, lat))[rowidxmin:rowidxmax,],
									raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt',
									radii = radii,
									fp = '/mnt/research/nasabio/temp',
									filetag = paste0('fia',slice))

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats/30m/stats_',slice,'.r'))