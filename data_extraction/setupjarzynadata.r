# Set up data for Jarzyna's occupancy model
# Group bbs data into 5 segments per route (10 stops each)
# Make list of arrays, 1 per year. Dimensions are n routes x 5 x n species
# Also get average elevation for each segment as a covariate.

load('/mnt/research/nasabio/data/bbs/bbsworkspace.r')

library(dplyr)

bbsbysegment <- bbscov %>%
  mutate(stopno = as.numeric(unlist(regmatches(Stop, gregexpr('[0-9]+',Stop)))),
         segment = ceiling(stopno/10)) %>%
  cbind(fixedbbsmat) 

names(bbsbysegment)[8:ncol(bbsbysegment)] <- paste0('sp', 1:(ncol(fixedbbsmat)))

bbsbysegment <- bbsbysegment %>%
  group_by(year, rteNo, segment) %>%
  summarize_at(8:ncol(bbsbysegment), .funs = function(x) any(x > 0))

save(bbsbysegment, file = '/mnt/research/nasabio/data/bbs/occmod_bysegment.r')

# Load bbs coordinates to get the elevations.
bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')
elevfile <- '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt'
elevfile_ak <- '/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/AlaskaSE_30m_dem.vrt'
library(raster)
elevraster <- raster(elevfile)
bbselevs <- extract(elevraster, with(bbsll, cbind(lon, lat)), method = 'simple')
elevraster_ak <- raster(elevfile_ak) # Might only be southeastern Alaska.
bbselevs_ak <- extract(elevraster_ak, with(bbsll, cbind(lon, lat)), method = 'simple')

# Get elevations from raster's alt data to use for the northern Alaskan routes.
# Get the entire USA altitude mask to use for northern AK.
usa_alt <- getData('alt', path='/mnt/research/plz-lab/NEON/external_data/raw_external_data/srtm', country='USA')
alaska_alt <- usa_alt[[2]]
bbselevs_ak2 <- extract(alaska_alt, with(bbsll, cbind(lon, lat)), method = 'simple')

# Combine the three elevation rasters.
bbselevs[is.na(bbselevs)] <- bbselevs_ak[is.na(bbselevs)]
bbselevs[is.na(bbselevs)] <- bbselevs_ak2[is.na(bbselevs)]

bbsll_elev <- cbind(bbsll, bbselevs)
save(bbsll_elev, file = '/mnt/research/nasabio/data/bbs/occmod_elevs.r')