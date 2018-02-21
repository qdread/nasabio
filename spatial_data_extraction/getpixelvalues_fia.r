# Get pixel values for FIA
# Done in parallel

source('/mnt/research/nasabio/code/loadfiaall.r')
coords <- fiacoords
slice <- as.numeric(Sys.getenv('PBS_ARRAYID'))
n_slices <- 150

# Get row indexes for the slice of coordinate matrix to be extracted.
rowidx <- round(seq(0,nrow(fiacoords),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]
coords <- coords[rowidxmin:rowidxmax, ]

# Get table of raster file locations
vartable <- read.csv('/mnt/research/nasabio/data/geodiv_table_for_pixel.csv', stringsAsFactors = FALSE)

scratch_path <- '/mnt/ls15/scratch/groups/nasabio/VRTs'

values <- list()

for (i in 1:(nrow(vartable))) {
  print(vartable$Dataset.name[i])
  values[[i]] <- matrix(NA, nrow = nrow(coords), ncol = vartable$N.layers[i])
  for (j1 in 1:nrow(coords)) {
    values[[i]][j1, ] <- system2('gdallocationinfo', args = paste('-wgs84 -valonly', file.path(scratch_path, vartable$File.name[i]), coords$lon[j1], coords$lat[j1]), stdout = TRUE)
    print(j1)
  }
}

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

names(values) <- correct_var_names
values <- cbind(PLT_CN = coords$PLT_CN, values)

write.csv(values, paste0('/mnt/research/nasabio/data/fia/allgeodiv_v2/fia_geo_by_point_', slice, '.csv'), row.names = FALSE)