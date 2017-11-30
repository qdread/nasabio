# Old and new elevation extraction. Debug extractBox()

# Extract just a single point from the raster.
# Check the old raster vs the "big" raster

library(raster)
x <- raster('conus_30m_dem_big.vrt')
x_small <- raster('conus_30m_dem.vrt')

fiapnw <- read.csv('~/data/pnw.csv')
testcoords <- fiapnw[10001:10010,]

# The below are identical so it doesn't seem the _big raster is the problem.
extract(x = x, y = testcoords[,c('ACTUAL_LON', 'ACTUAL_LAT')])
extract(x = x_small, y = testcoords[,c('ACTUAL_LON', 'ACTUAL_LAT')])

source('/mnt/research/nasabio/code/extractbox.r')

# Load precalculated distances
load('/mnt/research/nasabio/data/precalcdist/distlogical_fia_elev.r')

extractBox(coords = with(testcoords, cbind(ACTUAL_LON, ACTUAL_LAT))[1,,drop=FALSE],
           raster_file = 'conus_30m_dem_big.vrt',
           radius = 300,
           fp = '/mnt/research/nasabio/temp',
           filetags = 'fiatest1')

radii <- c(5, 10, 20, 30, 40, 50) # in km

stats_test <- statsByRadius('/mnt/research/nasabio/temp/bbox_fiatest1.tif', radii = radii, is_brick = FALSE, trig = FALSE, use_abs = FALSE)


# Comparison of old and new stats for same point
n_slices<-5000
slice<-990
rowidx <- round(seq(0,22532,length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]
# 4458

# Find lowest elev sd for new data
minrow <- which.min(edgeo$sd)
cnlow <- edgeo$CN[minrow]
which(fiacoords$PLT_CN == cnlow)
