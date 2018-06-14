# Diagnose issues with the new unfuzzed data. Which one(s) are fucked up?????

library(dplyr)
library(sp)

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

# Means
emean5 <- edgeo %>% filter(radius == 5) %>% arrange(mean)

plot(pnwmap)
with(emean5, points(x = ACTUAL_LON, y = ACTUAL_LAT, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(mean, breaks = 20))]))


# Try old one
olded <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_geodiversity_stats.csv')

datr5old <- olded %>% filter(radius == 5, variable == 'elevation') %>% arrange(sd)
datr100old <- olded %>% filter(radius == 100, variable == 'elevation') %>% arrange(sd)

plot(pnwmap)
with(datr5old, points(x = lon, y = lat, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(sd, breaks = 20))]))

plot(pnwmap)
with(datr100old, points(x = lon, y = lat, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(sd, breaks = 20))]))

emeanold <- olded %>% filter(radius == 5, variable == 'elevation') %>% arrange(mean)
plot(pnwmap)
with(emeanold, points(x = lon, y = lat, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(mean, breaks = 20))]))

# Plot with ggplot
library(cowplot)
edmapdat <- edgeo %>% filter(radius==rad, !is.na(sd)) %>% arrange(sd) %>% rename(lon = ACTUAL_LON, lat = ACTUAL_LAT)
edmapdat_OLD <- bdold %>% filter(radius == rad, !is.na(sd_elev)) %>% arrange(sd_elev)
fiamap_ed <- ggplot(edmapdat, 
                    aes(x = lon, y = lat)) +
  borders('world', 'canada', fill = 'gray70') +
  borders('world', 'usa', fill = 'gray70') +
  borders('state') +
  geom_point(aes(color = sd), size = 0.75) +
  geom_rect(data=westcoast_scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=westcoast_scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscaleelev + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')) +
  geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  ggtitle(paste(rad, 'km radius'))

fiamap_ed_OLD <- ggplot(edmapdat_OLD, 
                    aes(x = lon, y = lat)) +
  borders('world', 'canada', fill = 'gray70') +
  borders('world', 'usa', fill = 'gray70') +
  borders('state') +
  geom_point(aes(color = sd_elev), size = 0.75) +
  geom_rect(data=westcoast_scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=westcoast_scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscaleelev + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')) +
  geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  ggtitle(paste(rad, 'km radius'))

# Also look at BBS geodiv stats to make sure that they are not wrong too.
edbbs <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_geodiversity_stats.csv', stringsAsFactors = FALSE)

# Beta-diversity ----------------------------------------------------------

# Look at beta-diversity.
bd <- read.csv(file.path(fp, 'fia_betapartfinal_to100_wide.csv'))

bdgeo <- bd %>%
  mutate(radius = radius / 1000) %>%
  rename(CN = PLT_CN) %>%
  left_join(fiacoords) 


bdatr5 <- bdgeo %>% filter(radius == 5) %>% arrange(sorensen_taxonomic_total)
bdatr100 <- bdgeo %>% filter(radius == 100) %>% arrange(sorensen_taxonomic_total)

plot(pnwmap)
with(bdatr5, points(x = ACTUAL_LON, y = ACTUAL_LAT, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(bray_taxonomic_total, breaks = 20))]))

plot(pnwmap)
with(bdatr100, points(x = ACTUAL_LON, y = ACTUAL_LAT, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(bray_taxonomic_total, breaks = 20))]))

# Gamma-diversity ---------------------------------------------------------

gd <- read.csv(file.path(fp, 'fia_gammadiv.csv'))

gdgeo <- gd %>%
  rename(CN = PLT_CN) %>%
  left_join(fiacoords) 


gdatr5 <- gdgeo %>% filter(radius == 5) %>% arrange(shannon)
gdatr100 <- gdgeo %>% filter(radius == 100) %>% arrange(shannon)

plot(pnwmap)
with(gdatr5, points(x = ACTUAL_LON, y = ACTUAL_LAT, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(shannon, breaks = 20))]))

plot(pnwmap)
with(gdatr100, points(x = ACTUAL_LON, y = ACTUAL_LAT, pch = 19, cex = 0.5, col = my_colors(20)[as.numeric(cut(shannon, breaks = 20))]))

# Look at old gamma richness
gdold <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_gammadiv.csv')

# Join old and new gamma richness. Should be the same.