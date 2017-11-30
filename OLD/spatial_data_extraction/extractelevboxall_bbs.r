n_slices <- 2000 
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
namepart1 <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_'
namepart2 <- '_big.vrt'
varname <- c('dem','aspect','slope','TPI')
scratch_path <- Sys.getenv('TMPDIR')

# BBS lat long coordinates
bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')

# Load precalculated distances
load('/mnt/research/nasabio/data/precalcdist/distlogical_bbs_elev.r')

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(bbsll),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

stats_by_point <- list()
radii <- c(50, 75, 100, 150, 200, 300)

pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)
		   
for (j in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, j)
	stats_j <- list()
	for (k in 1:length(varname)) {
	istrig <- ifelse(varname[k] == 'aspect', TRUE, FALSE)
	isabs <- ifelse(varname[k] == 'TPI', TRUE, FALSE)
	extractBox(coords = with(bbsll, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = paste0(namepart1, varname[k], namepart2),
		   radius = 300,
		   fp = scratch_path,
		   filetags = paste('bbselev', row.names(bbsll)[j], k, sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_bbselev_', row.names(bbsll)[j], '_', k, '.tif')
	stats_j[[k]] <- statsByRadius(file_j, radii = radii, is_brick = FALSE, trig = istrig, use_abs = isabs)
	if (file.exists(file_j)) deleteBox(file_j)
}	
stats_by_point[[length(stats_by_point) + 1]] <- stats_j
}

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/bbs/elevstats/big30m/stats_',slice,'.r'))
