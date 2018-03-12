# Map of bird conservation region boundaries, state boundaries, and BBS route locations that we used, all on the same map

# Load BCR boundaries
library(sp)
library(rgdal)
library(ggplot2)
library(dplyr)
library(purrr)
fpregion <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/regions'
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
states <- read.csv('~/R/states_albers.csv', stringsAsFactors = FALSE)
bbscoords <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_correct_route_centroids.csv')

aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
bcr <- spTransform(bcr, CRSobj = CRS(aea_crs))

bcr@data <- bcr@data %>%
  mutate(id = rownames(bcr@data), region = as.character(BCRNAME))

region_fort <- bcr %>%
  subset(COUNTRY == 'USA' & !PROVINCE_S %in% c('ALASKA','HAWAII','HAWAIIAN ISLANDS')) %>%
  fortify(bcr, region = 'id') %>% 
  left_join(bcr@data, by = 'id')

blktheme <- theme_bw() + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'black'), 
        panel.border = element_blank(), 
        plot.background = element_rect(fill = 'black'), 
        legend.position = c(0.13,0.1), 
        legend.direction = 'horizontal', 
        legend.title = element_blank())

p <- ggplot(region_fort) +
    geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
    geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'gray50', size = 0.75) +
    geom_point(data = bbscoords %>% filter(lat < 50), aes(x=lon.1, y=lat.1), color = 'indianred', size = 0.5) +
    coord_equal() +
    blktheme
ggsave('C:/Users/Q/Dropbox/presentations/sesync2018/bbsbcrmap.png', p, height = 6, width = 9, dpi = 400)
