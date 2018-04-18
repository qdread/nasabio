# Get pixel values for BBS
# Done in parallel
# 14 Feb

# Edited 17 Apr: change from centroids to midpoints.

coords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv', stringsAsFactors=FALSE)

slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
n_slices <- 25

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(coords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]
coords <- coords[rowidxmin:rowidxmax, ]

# Get table of raster file locations
vartable <- read.csv('/mnt/research/nasabio/data/geodiv_table_for_pixel.csv', stringsAsFactors = FALSE)

scratch_path <- '/mnt/ls15/scratch/groups/nasabio/VRTs'


######################################
# Serial version

values <- list()

for (i in 1:(nrow(vartable))) {
  print(vartable$Dataset.name[i])
  values[[i]] <- matrix(NA, nrow = nrow(coords), ncol = vartable$N.layers[i])
  pb <- txtProgressBar(0, nrow(coords), style=3)
  for (j1 in 1:nrow(coords)) {
    setTxtProgressBar(pb, j1)
    values[[i]][j1, ] <- system2('gdallocationinfo', args = paste('-wgs84 -valonly', file.path(scratch_path, vartable$File.name[i]), coords$lon[j1], coords$lat[j1]), stdout = TRUE)
    
  }
  close(pb)
}

######################################

# Assemble output.
values <- do.call('cbind', values)
class(values) <- 'numeric'
values[values == -9999] <- NA
values[abs(values) > 1e100] <- NA
values <- as.data.frame(values)

# Give meaningful column names to values.
correct_var_names <- list('elevation_30m', 'elevation_30m_tri', 'elevation_30m_roughness',
						  'slope_30m', 
						  'aspect_sin_30m', 'aspect_cos_30m', 
						  paste0('bio', 1:19,'_1k'),
						  paste0('bio', 1:19,'_1k_tri'),
						  paste0('bio', 1:19,'_1k_roughness'),
						  paste0('biocloud', 1:8,'_1k'),
						  paste0('biocloud', 1:8,'_1k_tri'),
						  paste0('biocloud', 1:8,'_1k_roughness'),
						  c('dhi_fpar_1k', 'dhi_gpp_1k', 'dhi_lai8_1k', 'dhi_ndvi_1k'),
						  c('dhi_fpar_1k_tri', 'dhi_gpp_1k_tri', 'dhi_lai8_1k_tri', 'dhi_ndvi_1k_tri'),
						  c('dhi_fpar_1k_roughness', 'dhi_gpp_1k_roughness', 'dhi_lai8_1k_roughness', 'dhi_ndvi_1k_roughness'),
						  'human_footprint_1k', 'human_footprint_1k_tri', 'human_footprint_1k_roughness',
						  'nightlight_500m', 'nighlight_500m_tri', 'nightlight_500m_roughness',
						  'geological_age_1k', 'soil_type_5k')

names(values) <- unlist(correct_var_names)
values <- cbind(rteNo = coords$rteNo, values)

write.csv(values, paste0('/mnt/research/nasabio/data/bbs/allgeodiv_v2/bbs_geo_by_point_',slice,'.csv'), row.names = FALSE)


####################################

# Combine output.
bbs_geo <- lapply(1:25, function(slice) read.csv(paste0('/mnt/research/nasabio/data/bbs/allgeodiv_v2/bbs_geo_by_point_', slice, '.csv'), stringsAsFactors = FALSE))
bbs_geo <- lapply(bbs_geo, setNames, nm = c('rteNo', unlist(correct_var_names)))
bbs_geo <- do.call(rbind, bbs_geo)
write.csv(bbs_geo, '/mnt/research/nasabio/data/bbs/bbs_geo_by_point.csv', row.names = FALSE)