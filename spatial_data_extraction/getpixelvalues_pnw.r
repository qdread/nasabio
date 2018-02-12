# Get pixel values for FIA
# Done in parallel
# Edit 09 Feb for PNW

source('/mnt/research/nasabio/code/loadfiaall.r')

coords <- read.csv('~/data/pnw.csv', stringsAsFactors=FALSE)
names(coords)[5:6] <- c('lat','lon')

# Get table of raster file locations
vartable <- read.csv('/mnt/research/nasabio/data/geodiv_table_for_pixel.csv', stringsAsFactors = FALSE)

scratch_path <- file.path(Sys.getenv('SCRATCH'), 'geo')

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

# Assemble output.
values <- do.call('cbind', values)
class(values) <- 'numeric'
values[values == -9999] <- NA
values[abs(values) > 1e100] <- NA
values <- as.data.frame(values)

# Give meaningful column names to values.
var_names <- c(paste0('bio', 1:19,'_1k'), paste0('biocloud', 1:8, '_1k'), 
               'dhi_fpar_1k', 'dhi_gpp_1k', 'dhi_lai8_1k', 'dhi_ndvi_1k',
               'elevation_30m', 'slope_30m', 'human_footprint_1k', 'nightlight_500m',
               'geological_age_1k', 'soil_type_5k',
               paste0('TRI_bio', 1:19,'_1k'), paste0('TRI_biocloud', 1:8,'_1k'), 
               'TRI_elevation_30m', 
               'TRI_dhi_fpar_1k', 'TRI_dhi_gpp_1k', 'TRI_dhi_lai8_1k', 'TRI_dhi_ndvi_1k',
               'TRI_human_footprint_1k',
               'TRI_nightlight_500m')

names(values) <- var_names
values <- cbind(PLT_CN = coords$CN, values)

write.csv(values, '/mnt/research/nasabio/data/fia/fia_geo_by_point_correctedpnw.csv', row.names = FALSE)