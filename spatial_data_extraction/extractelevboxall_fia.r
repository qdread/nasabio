# Edit 10 Nov: Just do every plot separately!
n_slices <- 22531
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
namepart1 <- '/mnt/research/nasabio/data/dem/SRTM_30m_DEM/VRTs/conus_30m_'
namepart2 <- '_big.vrt'
varname <- c('dem','aspect','slope','TPI')
scratch_path <- Sys.getenv('TMPDIR')

# FIA lat long coordinates
source('/mnt/research/nasabio/code/loadfia.r')

# Function to do the extracting
source('/mnt/research/nasabio/code/extractbox.r')

# Load precalculated distances
load('/mnt/research/nasabio/data/precalcdist/distlogical_fia_elev.r')

# Get row indexes for the slice of coordinate matrix to be extracted.
j <- slice

stats_by_point <- list()
radii <- c(5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300) # in km

	stats_j <- list()
	for (k in 1:length(varname)) {
	istrig <- ifelse(varname[k] == 'aspect', TRUE, FALSE)
	isabs <- ifelse(varname[k] == 'TPI', TRUE, FALSE)
	extractBox(coords = with(fiacoords, cbind(lon, lat))[j,,drop=FALSE],
		   raster_file = paste0(namepart1, varname[k], namepart2),
		   radius = 300,
		   fp = scratch_path,
		   filetags = paste('fiaelev', row.names(fiacoords)[j], k, sep = '_'))
	file_j <- paste0(scratch_path, '/bbox_fiaelev_', row.names(fiacoords)[j], '_', k, '.tif')
	stats_j[[k]] <- statsByRadius(file_j, radii = radii, is_brick = FALSE, trig = istrig, use_abs = isabs)
	if (file.exists(file_j)) deleteBox(file_j)
	}
	stats_by_point[[length(stats_by_point) + 1]] <- stats_j

save(stats_by_point, file = paste0('/mnt/research/nasabio/data/fia/elevstats/unfuzzed/stats_',slice,'.r'))
