library(foreach)
library(doParallel)

registerDoParallel(cores = 25)


stats_by_point <- foreach(i = rowidxmin:rowidxmax) %dopar% {
	stats_i <- emptydf
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
				stats_i[r, 2:6] <- c(mean(eir), sd(eir), min(eir), max(eir), length(eir)) # Calculate summary stats
			}
		}
	}
	return(stats_i)
}
