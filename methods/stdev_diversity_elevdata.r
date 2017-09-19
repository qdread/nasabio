# proof of concept: richness of elevations is very similar to standard deviation of elevations

# Use FIA points. 100 points x radii
# Extract the radii with "easy" raster method
# Find standard deviation, as well as richness of the elevations rounded to the nearest whole number.

raster_file_name <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_dem.vrt'
scratch_path <- Sys.getenv('TMPDIR')


# FIA lat long coordinates
library(dplyr)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1]) %>%
  ungroup

fiacoordssample <- fiacoords[sample(10001:nrow(fiacoords), 100), ] %>% dplyr::select(lon, lat)

source('/mnt/research/nasabio/code/extractbox.r')

radii <- seq(5, 100, by = 5)
truncation_points <- c(0.01, 0.1, 1, 10, 100)

# Create precalculated distances.

lon1 <- -98.59595
lat1 <- 39.79797

extractBox(coords = cbind(lon1, lat1),
		   raster_file = elevfile,
		   radius = 100,
		   fp = '/mnt/research/nasabio/temp',
		   filetags = 'elev100')

# Get proximity map internally
	x <- raster('/mnt/research/nasabio/temp/bbox_elev100.tif')
	xvals <- as.data.frame(x)
	xcoords <- xyFromCell(x, 1:ncell(x))
	xdist <- spDistsN1(xcoords, c(lon1,lat1), longlat=TRUE)
	
	idx_r <- list()
	
	for (r in 1:length(radii)) {
		idx_r[[r]] <- xdist <= radii[r]
	}

	# Create logical matrix.
	# Use column r to subset xvalues each time.
	distr <- do.call('cbind', idx_r)



summ_stats <- function(boxfile, distr, radii, truncations) {
	require(raster)
	require(vegan)
	if (!file.exists(boxfile)) return(NA)
	x <- raster(boxfile)
	dat <- as.data.frame(x)
	res <- list()
	
	for (r in 1:ncol(distr)) {
		dat_r <- na.omit(dat[distr[, r], ])
		richnesses <- sapply(truncations, function(x) length(unique(round_any(dat_r, x))))
		diversities <- sapply(truncations, function(x) diversity(round_any(dat_r - min(dat_r), x)))
		sds <- sapply(truncations, function(x) sd(round_any(dat_r, x)))
		res[[r]] <- data.frame(radius = radii[r], truncation = truncations, richness = richnesses, diversity = diversities, sd = sds)
	}
	do.call('rbind', res)
}

stats_by_point <- list()

pb <- txtProgressBar(1, nrow(fiacoordssample), style=3)
		   
for (j in 1:nrow(fiacoordssample)) {
	setTxtProgressBar(pb, j)
	extractBox(coords = with(fiacoordssample, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = raster_file_name,
		   radius = 300,
		   fp = scratch_path,
		   filetags = paste('fiatest', row.names(fiacoords)[j], sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_fiatest_', row.names(fiacoords)[j], '.tif')
	stats_by_point[[length(stats_by_point) + 1]] <- summ_stats(boxfile = file_j, distr = idx_r, radii = radii, truncations = truncation_points)
	if (file.exists(file_j)) deleteBox(file_j)
}	

close(pb)

save(stats_by_point, file = '/mnt/research/nasabio/data/teststats.r')




	