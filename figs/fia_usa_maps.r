# FIA maps formatted for paper.
# Forked version 09 Jan: Entire USA (elevation only).

# NOTE: All needed CSVs are on the hpcc but I have downloaded them locally because it is faster. 
# Change file path to the second file path to get files from hpcc.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

# Load coordinates
fiacoords <- read.csv('~/FIA/FIA10nov/allfia.csv') %>% setNames(nm = c('PLT_CN', 'lat', 'lon'))

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_usa_elev_only.csv')) %>% left_join(fiacoords)
ad <- read.csv(file.path(fp, 'fiausa_alpha.csv')) %>% left_join(fiacoords)
bd <- read.csv(file.path(fp, 'fiausa_betatd.csv')) %>% left_join(fiacoords)
gd <- read.csv(file.path(fp, 'fiausa_gamma.csv')) %>% left_join(fiacoords)

# Correct errors in ed
ed <- filter(ed, min > -200, max < 5000)

library(dplyr)
library(cowplot)

# Customized function to add scale bars and north arrows to maps.
# See http://stackoverflow.com/questions/39067838/parsimonious-way-to-add-north-arrow-and-scale-bar-to-ggmap

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9))
colscalealpha <- scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 2)(9))
colscalegamma <- scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 2)(9))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = RColorBrewer::brewer.pal(9, 'YlOrRd'))


whitemaptheme <- theme_bw() + theme(panel.grid = element_blank(), 
                                      legend.position = c(0.3, 0.03),
                                      legend.background = element_rect(fill='transparent', color='transparent'),
                                      legend.direction = 'horizontal',
                                      axis.title = element_blank())

longn <- seq(-120, -70, by = 10)
latn <- seq(25, 50, by = 5)
fia_xs <- scale_x_continuous(breaks = longn, labels = paste0(abs(longn), '° W'))
fia_ys <- scale_y_continuous(breaks = latn, labels = paste0(latn, '° N'))

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps'

latbds = c(25,50)
lonbds <- c(-125, -67)

radii <- c(5, 10, 20, 50, 100)

single_fia_map <- function(dat, zvar, colscale, title, subtitle) {
  ggplot(dat, aes(x = lon, y = lat)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'mexico', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(aes_string(color = zvar), size = 0.1) + 
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    whitemaptheme + 
    theme(legend.position = c(0.1, 0.05)) +
    colscale + 
    panel_border(colour = 'black') +
    fia_xs + fia_ys +
    ggtitle(title, subtitle)
}

for (rad in radii) {
  
  bdmapdat <- bd %>% dplyr::filter(radius == rad, !is.na(beta_td_sorensen)) %>% arrange(beta_td_sorensen)
  admapdat <- ad %>% dplyr::filter(radius == rad, !is.na(shannon)) %>% mutate(shannon = exp(shannon)) %>% arrange(shannon)
  edmapdat <- ed %>% dplyr::filter(radius == rad, !is.na(sd), variable == 'elevation') %>% arrange(sd)
  gdmapdat <- gd %>% dplyr::filter(radius == rad, !is.na(shannon)) %>% mutate(shannon = exp(shannon)) %>% arrange(shannon)
  
  ad_map <- single_fia_map(dat = admapdat, zvar = 'shannon', colscale = colscalealpha, title = 'FIA taxonomic alpha-diversity', subtitle = paste('True Shannon diversity,', rad, 'km radius'))
  bd_map <- single_fia_map(dat = bdmapdat, zvar = 'beta_td_sorensen', colscale = colscalebeta, title = 'FIA taxonomic beta-diversity', subtitle = paste('Sorensen, abundance,', rad, 'km radius'))
  gd_map <- single_fia_map(dat = gdmapdat, zvar = 'shannon', colscale = colscalegamma, title = 'FIA taxonomic gamma-diversity', subtitle = paste('True Shannon diversity,', rad, 'km radius'))
  ed_map <- single_fia_map(dat = edmapdat, zvar = 'sd', colscale = colscaleelev, title = 'FIA plots elevation standard deviation', subtitle = paste(rad, 'km radius'))
  
  
  
  ggsave(file.path(fpfig, paste0('fia_usa_alpha_shannon_', rad, '_km.png')), ad_map, height = 8, width = 12, dpi = 400)
  ggsave(file.path(fpfig, paste0('fia_usa_beta_sorensen_', rad, '_km.png')), bd_map, height = 8, width = 12, dpi = 400)
  ggsave(file.path(fpfig, paste0('fia_usa_gamma_shannon_', rad, '_km.png')), gd_map, height = 8, width = 12, dpi = 400)
  ggsave(file.path(fpfig, paste0('fia_usa_stdev_elev_', rad, '_km.png')), ed_map, height = 8, width = 12, dpi = 400)
}

# Add the incidence-based metrics to maps too.

for (rad in radii) {
  
  bdmapdat <- bd %>% dplyr::filter(radius == rad, !is.na(beta_td_sorensen_pa)) %>% arrange(beta_td_sorensen_pa)
  admapdat <- ad %>% dplyr::filter(radius == rad, !is.na(richness)) %>% arrange(richness)
  gdmapdat <- gd %>% dplyr::filter(radius == rad, !is.na(richness)) %>% arrange(richness)
  
  ad_map <- single_fia_map(dat = admapdat, zvar = 'richness', colscale = colscalealpha, title = 'FIA taxonomic alpha-diversity', subtitle = paste('Richness,', rad, 'km radius'))
  bd_map <- single_fia_map(dat = bdmapdat, zvar = 'beta_td_sorensen_pa', colscale = colscalebeta, title = 'FIA taxonomic beta-diversity', subtitle = paste('Sorensen, incidence,', rad, 'km radius'))
  gd_map <- single_fia_map(dat = gdmapdat, zvar = 'richness', colscale = colscalegamma, title = 'FIA taxonomic gamma-diversity', subtitle = paste('Richness,', rad, 'km radius'))

  
  ggsave(file.path(fpfig, paste0('fia_usa_alpha_richness_', rad, '_km.png')), ad_map, height = 8, width = 12, dpi = 400)
  ggsave(file.path(fpfig, paste0('fia_usa_beta_sorensenincidence_', rad, '_km.png')), bd_map, height = 8, width = 12, dpi = 400)
  ggsave(file.path(fpfig, paste0('fia_usa_gamma_richness_', rad, '_km.png')), gd_map, height = 8, width = 12, dpi = 400)
}
