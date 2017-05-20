# Figures for FIA beta-diversity at different scales, regressed on elevational diversity

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/betadiv/'

bd <- read.csv(file.path(fp, 'fia_betatd.csv'), stringsAsFactors = FALSE)
ed <- read.csv(file.path(fp, 'fia_elev_stats_noalaska.csv'), stringsAsFactors = FALSE)

library(dplyr)

bd <- bd %>%
  mutate(radius = radius/1000) %>% # put radius in km
  left_join(ed)

# Scatterplots

library(cowplot)

fia_beta_plot <- ggplot(bd %>% subset(radius %in% c(5,10,20,50,100)), aes(x = sd_elev, y = beta_pairwise_abundance)) +
  geom_point(alpha = 0.33) +
  stat_smooth() +
  facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km')), scales = 'free_x') +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Pairwise abundance-weighted beta-diversity', expand = c(0,0)) +
  labs(x = 'Standard deviation of elevation')

ggsave(file.path(fpfig, 'fiabetadiversity.png'), fia_beta_plot, height = 9, width = 13.5, dpi = 400)

# Alpha diversity (taxonomic)

load(file.path(fp, 'fia_diversitymetrics.RData'))

# Also add functional diversity to this.
fia_fd <- read.csv(file.path(fp, 'fia_fd.csv'))
plot_diversity <- cbind(plot_diversity, fia_fd)

# Calculate plot diversity within radii

library(sp)

fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = F)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

plot_diversity <- left_join(plot_diversity, fiacoords)
radii <- c(5, 10, 20, 50, 100)


neighbordiv <- function(x) {
  neighbordists <- spDistsN1(pts = cbind(plot_diversity$lon, plot_diversity$lat), pt = c(x$lon, x$lat), longlat = TRUE)
  commdat <- list()
  for (i in 1:length(radii)) {
    neighbors_incircle <- plot_diversity[neighbordists <= radii[i], ]
    commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], richness = median(richness,na.rm=T), shannon_basalarea = median(shannon_basalarea, na.rm=T), mpd = median(mpd_z, na.rm=T), mntd = median(mntd_z, na.rm=T)))
  }
  as.data.frame(do.call('rbind', commdat))
}

fia_alpha <- plot_diversity %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>% 
  do(neighbordiv(.)) # Takes 4 minutes on my machine.


# Join with elevational diversity
ad <- left_join(fia_alpha, ed)

fia_alpha_plot <- ggplot(ad, aes(x = sd_elev, y = shannon_basalarea)) +
  geom_point(alpha = 0.33) +
  stat_smooth() +
  facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km')), scales = 'free_x') +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Median Shannon alpha-diversity', expand = c(0,0)) +
  labs(x = 'Standard deviation of elevation')

ggsave(file.path(fpfig, 'fiaalphadiversity.png'), fia_alpha_plot, height = 9, width = 13.5, dpi = 400)

####################################################

# Added 03 May: maps of elevational diversity, alpha-diversity, and beta-diversity at matching radii.

# First load bd and ad. Join both with ed.
# Do not use Alaska here.

latbds <- c(32.5, 50)
lonbds <- c(-125, -114)
aklatbds <- c(52, 61.5)
aklonbds <- c(-153.8, -125)
bothlatbds <- c(32.5, 61.5)
bothlonbds <- c(-153.8, -114)

purple9colors <- c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = purple9colors)	  
colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.3)(9))
colscalealpha <- scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.3)(9))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = RColorBrewer::brewer.pal(9, 'YlOrRd'))

blackmaptheme <- theme_void() + theme(panel.grid = element_blank(), 
                                      panel.background = element_rect(color = 'black', fill = 'black'), 
                                      plot.background = element_rect(color = 'black', fill = 'black'),  
                                      legend.position = c(0.25, 0.07), 
                                      legend.direction = 'horizontal',
                                      text = element_text(color = 'white'))




for (rad in radii) {
  bdmapdat <- bd %>% filter(radius == rad, !is.na(beta_pairwise_abundance)) %>% arrange(beta_pairwise_abundance)
  
  fiamap_bd <- ggplot(bdmapdat, 
                       aes(x = lon, y = lat, color = beta_pairwise_abundance)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(size = 0.75) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    blackmaptheme + colscalebeta +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('fia_map_',rad,'km_beta.png')
  ggsave(file.path(fpfig, fname), fiamap_bd, height = 8, width = 6, dpi = 400)
  
  admapdat <- ad %>% filter(radius == rad, !is.na(shannon_basalarea)) %>% arrange(shannon_basalarea)
  
  fiamap_ad <- ggplot(admapdat, 
                      aes(x = lon, y = lat, color = shannon_basalarea)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    #borders('state') +
    geom_point(size = 0.75) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    blackmaptheme + colscalealpha +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('fia_map_',rad,'km_alpha.png')
  ggsave(file.path(fpfig, fname), fiamap_ad, height = 8, width = 6, dpi = 400)
  
  edmapdat <- bd %>% filter(radius == rad, !is.na(sd_elev)) %>% arrange(sd_elev)
  
  fiamap_ed <- ggplot(edmapdat, 
                      aes(x = lon, y = lat, color = sd_elev)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    #borders('state') +
    geom_point(size = 0.75) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    blackmaptheme + colscaleelev +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('fia_map_',rad,'km_elevdiversity.png')
  ggsave(file.path(fpfig, fname), fiamap_ed, height = 8, width = 6, dpi = 400)
  
}




####################################################
# Added 03 May: BBS
# Edited 22 May: corrected metrics (tax div only)

bd <- read.csv(file.path(fp, 'bbs_beta_td_byroute.csv'), stringsAsFactors = FALSE)
ed <- read.csv(file.path(fp, 'bbs_elev_stats.csv'), stringsAsFactors = FALSE)

library(dplyr)

ed <- ed %>%
  #filter(year >= 2001, year <= 2011) %>%
  group_by(rteNo, radius) %>%
  summarize_all(.funs=mean, na.rm=T)

bd <- bd %>%
 # mutate(radius = radius/1000) %>% # put radius in km
  filter(year >= 2001, year <= 2011) %>%
  group_by(rteNo, radius) %>%
  summarize(beta_td_pairwise_presence=mean(beta_td_pairwise_pa, na.rm=T)) %>%
  left_join(ed)

library(cowplot)

radii <- c(50,75,100)

bbs_beta_plot <- ggplot(bd %>% filter(radius %in% radii), aes(x = sd_elev, y = beta_td_pairwise_presence)) +
  geom_point(alpha = 0.33) +
  stat_smooth() +
  facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km')), scales = 'free_x', nrow=2) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Pairwise abundance-weighted beta-diversity', expand = c(0,0)) +
  labs(x = 'Standard deviation of elevation')

ggsave(file.path(fpfig, 'bbsbetadiversity.png'), bbs_beta_plot, height = 9, width = 9, dpi = 400)

# Calculate alpha diversity within radii

bbs_alpha <- read.csv(file.path(fp,'bbs_alpha.csv')) # note: the lat and long in here are not good.


ad <- bbs_alpha %>%
  filter(year >= 2001, year <= 2011) %>%
  group_by(rteNo, radius) %>%
  summarize(richness=mean(richness, na.rm=T)) %>%
  left_join(ed)

bbs_alpha_plot <- ggplot(ad, aes(x = sd_elev, y = richness)) +
  geom_point(alpha = 0.33) +
  stat_smooth() +
  facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km')), scales = 'free_x', nrow=2) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Median richness', expand = c(0,0)) +
  labs(x = 'Standard deviation of elevation')
ggsave(file.path(fpfig, 'bbsalphadiversity.png'), bbs_alpha_plot, height = 9, width = 9, dpi = 400)


# 04 May: bbs maps.

latbds <- c(25, 50)
lonbds <- c(-125, -67)

#coords <- ad[,c('rteNo','radius','lon','lat')]
#bd <- left_join(bd,coords)

blackmaptheme <- theme_void() + theme(panel.grid = element_blank(), 
                                     panel.background = element_rect(color = 'black', fill = 'black'), 
                                     plot.background = element_rect(color = 'black', fill = 'black'),  
                                     legend.position = c(0.15, 0.07), 
                                     legend.direction = 'horizontal',
                                     text = element_text(color = 'white'))

radii <- c(50, 75, 100)
colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.8)(9))


for (rad in radii) {
  bbsmap_bd <- ggplot(bd %>% filter(radius == rad, !is.na(beta_td_pairwise_presence)) %>% arrange(beta_td_pairwise_presence), 
                      aes(x = lon, y = lat, color = beta_td_pairwise_presence)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'mexico', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(size = 1.5) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    blackmaptheme + colscalebeta +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('bbs_map_',rad,'km_beta.png')
  ggsave(file.path(fpfig, fname), bbsmap_bd, height = 6, width = 9, dpi = 400)

  bbsmap_ad <- ggplot(ad %>% filter(radius == rad, !is.na(richness)) %>% arrange(richness), 
                      aes(x = lon, y = lat, color = richness)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'mexico', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(size = 1.5) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    blackmaptheme + colscalealpha +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('bbs_map_',rad,'km_alpha.png')
  ggsave(file.path(fpfig, fname), bbsmap_ad, height = 6, width = 9, dpi = 400)
  
  bbsmap_ed <- ggplot(bd %>% filter(radius == rad, !is.na(sd_elev)) %>% arrange(sd_elev), 
                      aes(x = lon, y = lat, color = sd_elev)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'mexico', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(size = 1.5) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    blackmaptheme + colscaleelev +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('bbs_map_',rad,'km_elevdiversity.png')
  ggsave(file.path(fpfig, fname), bbsmap_ed, height = 6, width = 9, dpi = 400)
  
}
