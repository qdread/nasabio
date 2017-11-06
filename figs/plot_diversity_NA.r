# Diagnostic plot of new unfuzzed diversity to find locations of "NA" values

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_elev_stats_unfuzzed.csv'))
ad <- read.csv(file.path(fp, 'fia_alpha.csv'))
gd <- read.csv(file.path(fp, 'fia_gammadiv.csv'))

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_richness','beta_richness','gamma_richness')

library(dplyr)
library(sp)

biogeo <- ed %>%
  dplyr::select(PLT_CN, radius, sd) %>%
  filter(radius %in% radii) %>%
  rename(elevation_sd = sd) %>%
  left_join(ad %>% dplyr::select(PLT_CN, radius, richness) %>% rename(alpha_richness = richness)) %>%
  left_join(gd %>% dplyr::select(PLT_CN, radius, richness) %>% rename(gamma_richness = richness))

# Add latitude and longitudes from unfuzzed (on local drive only)
fiacoords <- read.csv('~/FIA/pnw.csv') %>% filter(complete.cases(.))

# Convert to albers
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
fia_aea <- spTransform(SpatialPoints(coords=fiacoords[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))


# Make some maps
biogeo <- cbind(biogeo, fia_aea@coords)

library(maps)
source('methods/spatial_fns.r')
pnwpoly <- make_map_polygons(states = c('California','Oregon','Washington'))
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
pnwpoly <- spTransform(pnwpoly, CRSobj = CRS(aea_crs))

plot(pnwpoly)
with(subset(biogeo, radius == 100), points(x = ACTUAL_LON, y = ACTUAL_LAT, col = c('blue', 'red')[is.na(gamma_richness) + 1]))
plot(pnwpoly)
with(subset(biogeo, radius == 100 & is.na(gamma_richness)), points(x = ACTUAL_LON, y = ACTUAL_LAT))
# They appear to be relatively random.

table(biogeo$radius, is.na(biogeo$elevation_sd)) # Same for each.

plot(pnwpoly)
with(subset(biogeo, radius == 100 & !is.na(elevation_sd)), points(x = ACTUAL_LON, y = ACTUAL_LAT))
# These also appear to be relatively random.