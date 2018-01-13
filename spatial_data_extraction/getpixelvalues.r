# Extract just pixel values for geo variables
# BBS and FIA
# QDR 12 Jan 2018
# Nasabioxgeo project

# Load coordinates
bbscoords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')

source('/mnt/research/nasabio/code/loadfiaall.r')

# Get table of raster file locations
vartable <- read.csv('/mnt/research/nasabio/data/geodiv_table_for_gdal.csv', stringsAsFactors = FALSE)

scratch_path <- file.path(Sys.getenv('SCRATCH'), 'geo')

bbsvalues <- list()

for (i in 1:(nrow(vartable)-2)) {
  print(vartable$Dataset.name[i])
  bbsvalues[[i]] <- matrix(NA, nrow = nrow(bbscoords), ncol = vartable$N.layers[i])
  for (j1 in 1:nrow(bbscoords)) {
    bbsvalues[[i]][j1, ] <- system2('gdallocationinfo', args = paste('-wgs84 -valonly', file.path(scratch_path, vartable$File.name[i]), bbscoords$lon[j1], bbscoords$lat[j1]), stdout = TRUE)
    print(j1)
  }
}

# Assemble output.
names(bbscoords)[4:5] <- c('lon_aea', 'lat_aea')
bbsvalues <- as.data.frame(cbind(bbscoords, do.call('cbind', bbsvalues)))

# Give meaningful column names to values.