# Set up data for Jarzyna's occupancy model
# Group bbs data into 5 segments per route (10 stops each)
# Make list of arrays, 1 per year. Dimensions are n routes x 5 x n species
# Also get average elevation for each segment as a covariate.

load('/mnt/research/nasabio/data/bbs/bbsworkspace.r') # This has already been pre-processed to remove nocturnal birds.

library(dplyr)

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

######
# Create objects for jags model.
load('/mnt/research/nasabio/data/bbs/occmod_elevs.r')

bbscov <- bbscov %>%
  mutate(rteNo = as.numeric(rteNo),
         stopno = as.numeric(unlist(regmatches(Stop, gregexpr('[0-9]+',Stop)))),
         segment = ceiling(stopno/10)) %>%
  left_join(bbsll_elev %>% select(rteNo, bbselevs) %>% rename(elev = bbselevs))

bbsjoin <- cbind(bbscov, fixedbbsmat)  
sp_names <- dimnames(fixedbbsmat)[[2]]

names(bbsjoin)[9:ncol(bbsjoin)] <- paste0('sp', 1:(ncol(fixedbbsmat)))

bbsbysegment <- bbsjoin %>%
  group_by(year, rteNo, elev, segment) %>%
  summarize_at(9:ncol(bbsjoin), .funs = function(x) any(x > 0))

#save(bbsbysegment, file = '/mnt/research/nasabio/data/bbs/occmod_bysegment.r')
#load('/mnt/research/nasabio/data/bbs/occmod_bysegment.r')

make_bbs_array <- function(df) {
  siteids <- unique(df$rteNo)
  res <- array(NA, dim = c(length(siteids), 5, ncol(df) - 3))
  for (site in 1:length(siteids)) {
    for (segment in 1:5) {
      res[site, segment, ] <- as.logical(df[df$rteNo == siteids[site] & df$segment == segment, -(1:3)])
    }
  }
  return(res)
}

bbs_arrays <- bbsbysegment %>%
  group_by(year) %>%
  do(x = make_bbs_array(.))

# By route, for initial values.
bbsbyroute <- bbsjoin %>%
  group_by(year, rteNo) %>%
  summarize_at(9:ncol(bbsjoin), .funs = function(x) any(x > 0))

df2mat <- function(x) as.matrix(x[,3:ncol(x)])

# Make this into a list.
bbsbyroutelist <- bbsbyroute %>%
  ungroup %>% group_by(year) %>%
  do(x = df2mat(.))

df2elevvector <- function(x) as.numeric(x$elev)

# Make elevations into a list as well.
bbselevvectors <- bbsjoin %>%
  group_by(year, rteNo) %>%
  do(x = df2elevvector(.))

#save(bbs_arrays, file = '/mnt/research/nasabio/data/bbs/occmod_arrays.r')

#load('/mnt/research/nasabio/data/bbs/occmod_bysegment.r')
#load('/mnt/research/nasabio/data/bbs/occmod_elevs.r')
#load('/mnt/research/nasabio/data/bbs/occmod_arrays.r')

# For the initial values we need the route observed presence/absences.
# Probably the same as the bbsmat byroute but I will recreate it here just to make sure that everything lines up.
#bbsbyroute <- cbind(bbscov, fixedbbsmat)

#names(bbsbyroute)[6:ncol(bbsbyroute)] <- paste0('sp', 1:(ncol(fixedbbsmat)))

bbsbyroute <- bbsbyroute %>%
  group_by(year, rteNo) %>%
  summarize_at(6:ncol(bbsbyroute), .funs = function(x) any(x > 0))

df2mat <- function(x) as.matrix(x[,3:ncol(x)])

# Make this into a list.
bbsbyroutelist <- bbsbyroute %>%
  ungroup %>% group_by(year) %>%
  do(x = df2mat(.))

# Match elevations with routes.
# We need a separate elevation vector for each year.
bbselevvectors <- bbscov %>% left_join(bbsll_elev %>% select(rteNo, bbselevs) %>% rename(elev = bbselevs)) %>%
  group_by(year) %>%
  do(x = function(.) as.numeric(.$elev))

# Save route ids and species ids to vectors, include in the final object so that no information gets mixed up.


