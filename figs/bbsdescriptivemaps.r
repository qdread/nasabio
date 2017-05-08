# Maps of BBS for powerpoint at NASA working group
# QDR 15 Feb 2017

# Map of BBS route locations, colored by different things.

fp <- '/mnt/research/aquaxterra/DATA/raw_data/BBS'
fp_fig <- '/mnt/research/nasabio/figs/descriptivemaps'
bd <- read.csv(file.path(fp, 'bbs_div_byroute.csv'), stringsAsFactors = FALSE)
bwgs <- read.csv(file.path(fp, 'bbs_wgs84_coords_byroute.csv'))

bd <- cbind(bd, bwgs)

library(dplyr)
library(ggplot2)
library(maptools)

# US Map with bbs routes.

dat <- bd %>%
	filter(year >= 2001, year <= 2011, !is.na(x_wgs84)) %>%
	group_by(rteNo) %>%
	summarize(lat = mean(y_wgs84),
			  lon = mean(x_wgs84),
			  richness = mean(richness))
			  

blackmaptheme <- theme_void() + theme(panel.grid = element_blank(), 
									  panel.background = element_rect(color = 'black', fill = 'black'), 
									  plot.background = element_rect(color = 'black', fill = 'black'),  
									  legend.position = c(0.2, 0.12), 
									  legend.direction = 'horizontal',
									  text = element_text(color = 'white'))
us_x <- scale_x_continuous(limits = c(-125, -67))
us_y <- scale_y_continuous(limits = c(25, 50))		

purple7colors <- c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')[-(1:2)]
colscale <- scale_colour_gradientn(name = 'Bird species richness', colours = purple7colors)	  
			  
bbsmap1 <- ggplot(dat) +
	borders('state', fill = 'gray90') +
	geom_point(aes(x = lon, y = lat, color = richness), size = 2) +
	coord_map(projection = 'albers', lat0 = 23, lat1 = 29.5) +
	us_x + us_y + colscale +
	blackmaptheme

ggsave(file.path(fp_fig, 'bbsroutes_richness.png'), bbsmap1, height = 5, width = 8, dpi = 400)

# Added 18 April: beta-diversity maps
betatd <- read.csv('/mnt/research/nasabio/data/bbs/bbs_beta_td.csv', stringsAsFactors=F)

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
latlong_crs <- '+proj=longlat +ellps=WGS72 +no_defs'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

bbs_aea <- SpatialPoints(coords = with(betatd, data.frame(x = x_aea, y = y_aea)), proj4string = CRS(aea_crs))
bbs_wgs <- spTransform(bbs_aea, CRSobj = CRS(wgs_crs))

betatd$lat <- bbs_wgs@coords[,2]
betatd$lon <- bbs_wgs@coords[,1]

dat <- betatd %>%
  filter(year >= 2001, year <= 2011, !is.na(lat)) %>%
  group_by(rteNo, Stop, radius) %>%
  summarize(x = mean(lon),
            y = mean(lat),
            beta_td = mean(beta_td_pairwise_presence))

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta diversity', colours = rainbow(10))	  

bbsbeta1000 <- ggplot(dat %>% filter(radius==1000)) +
  borders('state', fill = 'gray90') +
  geom_point(aes(x = x, y = y, color = beta_td), size = 2) +
  coord_map(projection = 'albers', lat0 = 23, lat1 = 29.5) +
  us_x + us_y + colscalebeta +
  blackmaptheme
bbsbeta2000 <- ggplot(dat %>% filter(radius==2000)) +
  borders('state', fill = 'gray90') +
  geom_point(aes(x = x, y = y, color = beta_td), size = 2) +
  coord_map(projection = 'albers', lat0 = 23, lat1 = 29.5) +
  us_x + us_y + colscalebeta +
  blackmaptheme
bbsbeta3000 <- ggplot(dat %>% filter(radius==3000)) +
  borders('state', fill = 'gray90') +
  geom_point(aes(x = x, y = y, color = beta_td), size = 2) +
  coord_map(projection = 'albers', lat0 = 23, lat1 = 29.5) +
  us_x + us_y + colscalebeta +
  blackmaptheme
bbsbeta4000 <- ggplot(dat %>% filter(radius==4000)) +
  borders('state', fill = 'gray90') +
  geom_point(aes(x = x, y = y, color = beta_td), size = 2) +
  coord_map(projection = 'albers', lat0 = 23, lat1 = 29.5) +
  us_x + us_y + colscalebeta +
  blackmaptheme
bbsbeta5000 <- ggplot(dat %>% filter(radius==5000)) +
  borders('state', fill = 'gray90') +
  geom_point(aes(x = x, y = y, color = beta_td), size = 2) +
  coord_map(projection = 'albers', lat0 = 23, lat1 = 29.5) +
  us_x + us_y + colscalebeta +
  blackmaptheme

ggsave(file.path(fp_fig, 'bbsbeta1km.png'), bbsbeta1000 + ggtitle('1 km radius'), height = 5, width = 8, dpi = 400)
ggsave(file.path(fp_fig, 'bbsbeta2km.png'), bbsbeta2000 + ggtitle('2 km radius'), height = 5, width = 8, dpi = 400)
ggsave(file.path(fp_fig, 'bbsbeta3km.png'), bbsbeta3000 + ggtitle('3 km radius'), height = 5, width = 8, dpi = 400)
ggsave(file.path(fp_fig, 'bbsbeta4km.png'), bbsbeta4000 + ggtitle('4 km radius'), height = 5, width = 8, dpi = 400)
ggsave(file.path(fp_fig, 'bbsbeta5km.png'), bbsbeta5000 + ggtitle('5 km radius'), height = 5, width = 8, dpi = 400)
