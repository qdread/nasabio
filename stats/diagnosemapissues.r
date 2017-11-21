# Diagnose issues with the new unfuzzed data. Which one(s) are fucked up?????

# Is it elevation?
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_elev_stats_unfuzzed.csv'))

# Join elevation with plot locations.

# Add latitude and longitudes from unfuzzed (on local drive only) 
# Also get rid of Alaska b/c no one cares and it isn't really a state.
fiacoords <- read.csv('~/FIA/pnw.csv') %>% filter(complete.cases(.))

# Project to Albers for drawing.
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
fia_aea <- spTransform(SpatialPoints(coords=fiacoords[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))

fiacoords <- transform(fiacoords, lon_aea = fia_aea@coords[,1], lat_aea = fia_aea@coords[,2])

edgeo <- ed %>%
  rename(CN = PLT_CN) %>%
  left_join(fiacoords) %>%
  filter(!is.na(sd))

# Draw maps. (quick version with base plot)
source('~/GitHub/nasabio/stats/spatial_fns.r')
pnwmap <- make_map_polygons(states = c('California','Oregon','Washington'))

my_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'))

datr5 <- edgeo %>% filter(radius == 5) %>% arrange(sd)
datr100 <- edgeo %>% filter(radius == 100) %>% arrange(sd)

plot(pnwmap)
with(datr5, points(x = ACTUAL_LON, y = ACTUAL_LAT, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(sd, breaks = 20))]))

plot(pnwmap)
with(datr100, points(x = ACTUAL_LON, y = ACTUAL_LAT, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(sd, breaks = 20))]))

# check this.
dat5100 <- datr5 %>% dplyr::select(CN, sd) %>% rename(sd5 = sd) %>% left_join(datr100 %>% dplyr::select(CN, sd) %>% rename(sd100 = sd)) 

with(dat5100, plot(sd5, sd100))

# Try old one
olded <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_geodiversity_stats.csv')

datr5old <- olded %>% filter(radius == 5, variable == 'elevation') %>% arrange(sd)
datr100old <- olded %>% filter(radius == 100, variable == 'elevation') %>% arrange(sd)

plot(pnwmap)
with(datr5old, points(x = lon, y = lat, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(sd, breaks = 20))]))

plot(pnwmap)
with(datr100old, points(x = lon, y = lat, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(sd, breaks = 20))]))
