# New: now that Kyla has pulled the entire USA in a single VRT file, all I need to do is get the maximally sized square box, get distances from each cell, then subset by radius for any radius we want.
# This version is for FIA.

n_slices <- 100 # Divide this into 100 little jobs, because this will take a while.
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))

# BBS lat long coordinates

library(dplyr)

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])


rad <- 1e5 # Maximal radius for bbs and fia, as anything above 100km has WAY too many pixels in it.

# Find row index of slice.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

library(raster)

# Load USA raster
usa <- raster('/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt')
usaext <- extent(usa)

# Radii of circles where we want summary stats
radii <- c(1, 5, 10, 20, 50) # in km. Just do up to 50 km for fia since there are more plots in that radius.

emptydf <- data.frame(radius = radii, mean_elev = NA, sd_elev = NA, min_elev = NA, max_elev = NA, n = NA)

# Data structure to hold the stats
stats_by_point <- replicate(length(rowidxmin:rowidxmax), emptydf, simplify = FALSE)

pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)

for (i in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb,i)
	if (!is.na(fiacoords$lat[i])) {
		loni <- fiacoords$lon[i]
		lati <- fiacoords$lat[i]
		if (loni > usaext@xmin & loni < usaext@xmax & lati > usaext@ymin & lati < usaext@ymax) {
			ei <- extract(usa, cbind(x=loni,y=lati), buffer=rad, cellnumbers = TRUE) # The actual extraction of elevations, will be slow
			ex <- xFromCell(usa, ei[[1]][,1]) # Find longitudes of the extracted elevations
			ey <- yFromCell(usa, ei[[1]][,1]) # Find latitudes of the extracted elevations
			disti <- spDistsN1(cbind(ex,ey), c(loni,lati), longlat=TRUE) # Distance from each cell to point i
			
			for (r in 1:length(radii)) {
				eir <- na.omit(ei[[1]][disti <= radii[r], 2]) # Elevations within the circle
				stats_by_point[[(i - rowidxmin + 1)]][r, 2:6] <- c(mean(eir), sd(eir), min(eir), max(eir), length(eir)) # Calculate summary stats
			}
		}
	}
	
}

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats/30m/stats_',slice,'.csv')