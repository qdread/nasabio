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


# Add bbs data to map -----------------------------------------------------

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

# BBS
bbsbio <- read.csv(file.path(fp, 'bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

bbsbio <- bbsbio %>%
  dplyr::select(rteNo, alpha_richness_100, beta_td_sorensen_pa_100, gamma_richness_100)
bbsgeo <- bbsgeo %>%
  dplyr::select(rteNo, elevation_5k_100_sd, bio1_5k_100_mean, geological_age_5k_100_diversity)

bbscoords <- bbscoords %>% left_join(bbsbio) %>% left_join(bbsgeo)
cscale <- function(b, n) scale_color_gradientn(name = n, colours = colorRampPalette(RColorBrewer::brewer.pal(9,'YlOrRd'), bias = b)(50))
biases <- c(1, 2, 0.8, 2, 1, 1)
legnames <- c('Alpha diversity', 'Beta diversity', 'Gamma diversity', 'Elevation variability', 'Temperature mean', 'Geological age diversity')

for (i in 6:11) {
dati <- bbscoords %>% 
  filter(lat < 50, !is.na(bbscoords[,i]))
p <- ggplot(region_fort) +
  geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'gray50', size = 0.75) +
  geom_point(data = dati, aes_string(x='lon.1', y='lat.1', color = names(bbscoords)[i]), size = 1) +
  coord_equal() +
  blktheme + cscale(biases[i-5], legnames[i-5]) + 
  theme(legend.text = element_text(color='white'),
        legend.title = element_text(color = 'white'),
        legend.background = element_rect(color = 'black', fill = 'black'))
ggsave(paste0('C:/Users/Q/Dropbox/presentations/sesync2018/bbs', names(bbscoords)[i], 'map.png'), p, height = 6, width = 9, dpi = 400)
}

# Map without BCRs with bigger dots
for (i in 6:8) {
  dati <- bbscoords %>% 
    filter(lat < 50, !is.na(bbscoords[,i]))
  p <- ggplot() +
    #geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
    geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'white', size = 0.75) +
    geom_point(data = dati, aes_string(x='lon.1', y='lat.1', color = names(bbscoords)[i]), size = 2.5) +
    coord_equal() +
    blktheme + cscale(biases[i-5], legnames[i-5]) + 
    theme(legend.text = element_text(color='white'),
          legend.title = element_text(color = 'white'),
          legend.background = element_rect(color = 'black', fill = 'black'))
  ggsave(paste0('C:/Users/Q/google_drive/NASABiodiversityWG/Conferences/IALE2018/bbs', names(bbscoords)[i], 'map.png'), p, height = 6, width = 9, dpi = 400)
}

# Same map for FIA --------------------------------------------------------

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
fiacoords <- read.csv('~/FIA/fiafuzzedbyQ.csv')

fiabio <- read.csv(file.path(fp, 'fia_unfuzzed/fia_bio_formixedmodels.csv'), stringsAsFactors = FALSE)
fiageo <- read.csv(file.path(fp, 'fia_unfuzzed/fia_geo_formixedmodels.csv'), stringsAsFactors = FALSE)

fiabio <- fiabio %>%
  dplyr::select(PLT_CN, alpha_richness, beta_td_sorensen_pa, gamma_richness)
fiageo <- fiageo %>%
  dplyr::select(PLT_CN, elevation_5k_100_sd, bio1_5k_100_mean, geological_age_5k_100_diversity)

fiacoords <- fiacoords %>% left_join(fiabio) #%>% left_join(fiageo)
cscale <- function(b, n) scale_color_gradientn(name = n, colours = colorRampPalette(RColorBrewer::brewer.pal(9,'YlOrRd'), bias = b)(50))
biases <- c(1, 0.8, 0.8, 2, 1, 1)
legnames <- c('Alpha diversity', 'Beta diversity', 'Gamma diversity', 'Elevation variability', 'Temperature mean', 'Geological age diversity')

for (i in 6:8) {
  dati <- fiacoords %>% 
    filter(latfuzz < 50, !is.na(fiacoords[,i]))
p <- ggplot() +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'white') +
  geom_point(data = dati, aes_string(x='lonfuzz_aea', y='latfuzz_aea', color = names(fiacoords)[i]), size = 1) +
  coord_equal() +
  cscale(biases[i-5], legnames[i-5]) +
  blktheme + 
  theme(legend.text = element_text(color='white'),
        legend.title = element_text(color = 'white'),
        legend.background = element_rect(color = 'black', fill = 'black'))
ggsave(paste0('C:/Users/Q/google_drive/NASABiodiversityWG/Conferences/IALE2018/fia', names(fiacoords)[i], 'map.png'), p, height = 6, width = 9, dpi = 400)
}