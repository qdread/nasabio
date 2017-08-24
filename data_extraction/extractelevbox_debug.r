# Test bounding box creation.

# FIA lat long coordinates
library(dplyr)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])
			
set.seed(517)
idx_test <- sample(nrow(fiacoords), 5)

extractBox(coords = with(fiacoords, cbind(lon, lat))[idx_test,],
		   raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt',
		   radius = 300,
		   fp = '/mnt/research/nasabio/temp',
		   filetags = paste('fiatest', row.names(fiacoords)[idx_test], sep = '_'))

stats_test <- list()
		   
for (i in 1:length(idx_test)) {
	print(i)
	file_i <- paste0('/mnt/research/nasabio/temp/bbox_fiatest_', row.names(fiacoords)[idx_test[i]], '.tif')
	stats_test[[i]] <- statsByRadius(file_i)
	if (file.exists(file_i)) deleteBox(file_i)
}	

	
### debug

coords = with(fiacoords, cbind(lon, lat))[idx_test,]
		   raster_file = '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt'
		   radius = 300
		   fp = '/mnt/research/nasabio/temp'
		   filetag = 'fiatest'
		   i=1

library(sp)
library(rgdal)

	latLongBox <- function(lon1, lat1, r) {
		# length of 1 degree latitude anywhere, and 1 degree longitude at equator, in km
		length1 <- 6371*2*pi/360
		latmin <- lat1 - r/length1
		latmax <- lat1 + r/length1
		lonmin <- lon1 - cos(lat1 * pi/180) * r/length1
		lonmax <- lon1 + cos(lat1 * pi/180) * r/length1
		cbind(c(lonmin, lonmin, lonmax, lonmax), c(latmin, latmax, latmax, latmin)) 
	}
	
	# define function to make the polygon into a json
	polygon2json <- function(p) {
		paste('{"type": "Polygon", "coordinates": [ [', paste(apply(p, 1, function(x){paste("[",x[1],", ",x[2],"]",sep="")}), collapse=","), '] ]}', sep="")
	}
	
loni <- coords[i,1]
lati <- coords[i,2]


box_r <- latLongBox(loni, lati, radius)
				box_json <- polygon2json(box_r)
				box_geojson <- readOGR(box_json, "OGRGeoJSON", verbose = FALSE,  p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
				# write the geojson to .shp
				writeOGR(box_geojson, fp, paste("temp_bbox", i, filetag, sep="_"), driver="ESRI Shapefile")
				
				
library(raster)
x <- raster(file.path(fp, 'bbox_1_fiatest.tif'))

# Get proximity map internally
xvals <- as.data.frame(x)
xcoords <- xyFromCell(x, 1:ncell(x))
disti <- spDistsN1(xcoords, c(loni,lati), longlat=TRUE)

radii <- c(1,5,10,20,50,100,150,200,300)
stats_i <- list()

for (r in 1:length(radii)) {
				print(radii[r])
				xvi <- na.omit(xvals[disti <= radii[r], ]) # Elevations within the circle
				stats_i[[r]] <- c(mean(xvi), sd(xvi), min(xvi), max(xvi), length(xvi)) # Calculate summary stats
			}

# Use precalculated proximity map with the following radii: c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300)

statsByRadius(file.path(fp, 'bbox_1_fiatest.tif'))


statsByRadius <- function(boxfile) {
	radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300)
	if (!file.exists(boxfile)) return(NA)
	x <- raster(boxfile)
	xvals <- as.data.frame(x)
	stats_r <- list()
	for (r in 1:ncol(idx_r)) {
		print(r)
		vals_r <- xvals[idx_r[, r], , drop = FALSE]
		means <- apply(vals_r, 2, mean, na.rm = TRUE)
		sds <- apply(vals_r, 2, sd, na.rm = TRUE)
		mins <- apply(vals_r, 2, min, na.rm = TRUE)
		maxes <- apply(vals_r, 2, max, na.rm = TRUE)
		nums <- apply(vals_r, 2, function(z) sum(!is.na(z)))
		stats_r[[r]] <- data.frame(radius = radii[r],
								   variable = names(xvals),
								   mean = means,
								   sd = sds,
								   min = mins,
								   max = maxes,
								   n = nums)
	}
	do.call('rbind', stats_r)	
}
