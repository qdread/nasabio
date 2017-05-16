# Some code to test the different extraction functions I developed.
# Must run on cluster due to the high ram demand.

library(dplyr)

source('~/code/fia/elev/extractfrombuffer_subsample.r')

# Load fia coordinates.
fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

# Radii of circles where we want summary stats
#radii <- c(0.01, 0.1, 1) 

#testcoords <- as.matrix(fiacoords[10001:10005,c('lon','lat')])

#extract_from_buffer(coords = testcoords, raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt', radii = radii)

#extract_from_buffer_subsample(coords = testcoords, raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt', radius = 5, samplesize = 10)

# Test how many subsamples is required to asymptotically reach a value.

testradius <- 500 # Large 500 km radius.
testcoords <- as.matrix(fiacoords[c(10001, 13295, 14002, 16454, 20100), c('lon', 'lat')])
testsamplesizes <- 10^(3:8)

summary_stats <- list()

for (i in 1:length(testsamplesizes)) {
	summary_stats[[i]] <- extract_from_buffer_subsample(coords = testcoords, 
														raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt', 
														radius = testradius, 
														samplesize = testsamplesizes[i])
}

save(summary_stats, file = '/mnt/research/nasabio/data/fia/testsummarystats.r')