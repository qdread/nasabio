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
