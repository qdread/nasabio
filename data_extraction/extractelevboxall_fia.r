n_slices <- 5000
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
namepart1 <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_'
namepart2 <- '_big.vrt'
varname <- c('dem','aspect','slope','TPI')
scratch_path <- Sys.getenv('TMPDIR')

# FIA lat long coordinates
library(raster)
source('/mnt/research/nasabio/code/loadfia.r')

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')

# Load precalculated distances
load('/mnt/research/nasabio/data/precalcdist/distlogical_fia_elev.r')

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

stats_by_point <- list()

pb <- txtProgressBar(rowidxmin, rowidxmax, style=3)
		   
for (j in rowidxmin:rowidxmax) {
	setTxtProgressBar(pb, j)
	stats_j <- list()
	for (k in 1:length(varname)) {
	istrig <- ifelse(varname == 'aspect', TRUE, FALSE)
	extractBox(coords = with(fiacoords, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = paste0(namepart1, varname[k], namepart2),
		   radius = 300,
		   fp = scratch_path,
		   filetags = paste('fiaelev', row.names(fiacoords)[j], k, sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_fiaelev_', row.names(fiacoords)[j], '_', k, '.tif')
	stats_j[[k]] <- statsByRadius(file_j, radii = radii, is_brick = FALSE, trig = istrig)
	if (file.exists(file_j)) deleteBox(file_j)
	}
	stats_by_point[[length(stats_by_point) + 1]] <- stats_j
}	

close(pb)

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats/big30m/stats_',slice,'.r'))
