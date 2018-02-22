# FIA maps for book chapter
# Same format as the BBS maps previously in the book chapter
# QDR 22 Feb 2018

fig_h <- 11 - 2*(2.5/2.54)
fig_w <- 8.5 - 2*(2.5/2.54)

# Must use eastern United States only.

library(dplyr)
library(cowplot)

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

# Load coordinates (TOP SECRET)
fiacoords <- read.csv('~/FIA/FIA10nov/allfia.csv') %>% setNames(nm = c('PLT_CN', 'lat', 'lon'))

# Exclude Western states
western_states <- c('California','Oregon','Washington', 'Idaho','Nevada','Arizona','New Mexico','Colorado','Utah','Wyoming','Montana')
source('stats/spatial_fns.r')
west_poly <- make_map_polygons(states = western_states)
fia_west <- SpatialPoints(fiacoords[,c('lon','lat')], proj4string = CRS(proj4string(west_poly))) %over% west_poly
fiacoords$east <- is.na(fia_west) 

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_usa_elev_only.csv')) %>% left_join(fiacoords)
ad <- read.csv(file.path(fp, 'fiausa_alpha.csv')) %>% left_join(fiacoords)
bd <- read.csv(file.path(fp, 'fiausa_betatd.csv')) %>% left_join(fiacoords)
gd <- read.csv(file.path(fp, 'fiausa_gamma.csv')) %>% left_join(fiacoords)

# Correct errors in ed
ed <- filter(ed, min > -200, max < 5000)

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9), breaks = c(0,0.5,1), labels = c(0,0.5,1))
colscalealpha <- scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 2)(9))
colscalegamma <- scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 2)(9))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 2)(9), breaks = c(0, 250, 500), labels = c(0, 250, 500))


whitemaptheme <- theme_bw() + theme(panel.grid = element_blank(), 
                                    legend.position = c(0.3, 0.03),
                                    legend.background = element_rect(fill='transparent', color='transparent'),
                                    legend.direction = 'horizontal',
                                    axis.title = element_blank())

longn <- seq(-105, -70, by = 10)
latn <- seq(25, 50, by = 5)
fia_xs <- scale_x_continuous(breaks = longn, labels = paste0(abs(longn), '° W'))
fia_ys <- scale_y_continuous(breaks = latn, labels = paste0(latn, '° N'))

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps/bookchapter'

latbds = c(25,50)
lonbds <- c(-105, -67)

rad <- 50 # Try 50 km.

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

# Make the data more granular so that there are no true coordinates included.

bdmapdat <- bd %>% 
  dplyr::filter(east, radius == rad, !is.na(beta_td_sorensen)) %>% 
  mutate(lat_round = plyr::round_any(lat, 0.1), lon_round = plyr::round_any(lon, 0.1)) %>%
  group_by(lat_round, lon_round) %>%
  summarize(beta_td_sorensen = mean(beta_td_sorensen)) %>%
  rename(lat = lat_round, lon = lon_round) %>%
  arrange(beta_td_sorensen)
admapdat <- ad %>% 
  dplyr::filter(east, radius == rad, !is.na(shannon)) %>% 
  mutate(shannon = exp(shannon), lat_round = plyr::round_any(lat, 0.1), lon_round = plyr::round_any(lon, 0.1)) %>% 
  group_by(lat_round, lon_round) %>%
  summarize(shannon = mean(shannon)) %>%
  rename(lat = lat_round, lon = lon_round) %>%
  arrange(shannon)
edmapdat <- ed %>% 
  dplyr::filter(east, radius == rad, !is.na(sd), variable == 'elevation', sd <= 1000) %>% 
  mutate(lat_round = plyr::round_any(lat, 0.1), lon_round = plyr::round_any(lon, 0.1)) %>%
  group_by(lat_round, lon_round) %>%
  summarize(sd = mean(sd)) %>%
  rename(lat = lat_round, lon = lon_round) %>%
  arrange(sd)
gdmapdat <- gd %>% dplyr::filter(east, radius == rad, !is.na(shannon)) %>% mutate(shannon = exp(shannon)) %>% arrange(shannon)

ad_map <- single_fia_map(dat = admapdat, zvar = 'shannon', colscale = colscalealpha, title = 'FIA taxonomic alpha-diversity', subtitle = paste('True Shannon diversity,', rad, 'km radius'))
bd_map <- single_fia_map(dat = bdmapdat, zvar = 'beta_td_sorensen', colscale = colscalebeta, title = 'FIA taxonomic beta-diversity', subtitle = paste('Sorensen, abundance,', rad, 'km radius'))
gd_map <- single_fia_map(dat = gdmapdat, zvar = 'shannon', colscale = colscalegamma, title = 'FIA taxonomic gamma-diversity', subtitle = paste('True Shannon diversity,', rad, 'km radius'))
ed_map <- single_fia_map(dat = edmapdat, zvar = 'sd', colscale = colscaleelev, title = 'FIA plots elevation standard deviation', subtitle = paste(rad, 'km radius'))


# Maps formatted for book -------------------------------------------------

source('figs/bbs_map_drawing_fns.r')

whitemaptheme <- theme_bw() + theme(panel.grid = element_blank(), 
                                    legend.position = c(0.1, 0.05),
                                    legend.background = element_rect(fill='transparent', color='transparent'),
                                    legend.direction = 'horizontal',
                                    axis.title = element_blank(),
                                    strip.background = element_blank())

deglab <- function(n, d) {
  paste0(as.character(abs(n)), '° ', d)
}



#### Edited map drawing fn
draw_fia_map_book <- function(dat, zvar, colscale, 
                              scalebar, 
                              latbds = c(25,50), lonbds = c(-100, -67), 
                              maptheme = whitemaptheme, 
                              xsc = scale_x_continuous(breaks = seq(-100, -70, by=10), labels = deglab( seq(-100, -70, by=10), 'W')), 
                              ysc = scale_y_continuous(breaks = c(25,35,45), labels = deglab(c(25, 35, 45), 'N')), 
                              the_arrow,
                              northlabel, 
                              maptitle = NULL,
                              fp = '~', fname = 'plot.png', write_to_file = TRUE, img_h = 7, img_w = 12) {
  the_map <- ggplot(dat, aes(x = lon, y = lat)) +
    borders('world', 'canada', fill = 'gray90') +
    borders('world', 'mexico', fill = 'gray90') +
    borders('world', 'usa', fill = 'gray90') +
    borders('state') +
    geom_point(aes_string(color = zvar), size = 0.25) +
    geom_rect(data=scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill.col), inherit.aes=F,
              show.legend = F,  color = "black") +
    scale_fill_manual(values = c('black', 'white')) +
    geom_text(data=scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 1) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    maptheme + colscale + panel_border(colour = 'black') +
    xsc + ysc + 
    the_arrow + 
    northlabel + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          legend.position = 'bottom')
  if (!is.null(maptitle)) the_map <- the_map + ggtitle(maptitle)
  if (write_to_file) ggsave(file.path(fp, fname), the_map, height = img_h, width = img_w, dpi = 400) 
  if (!write_to_file) return(the_map)
}

small_low_arrow <- geom_line(data = data.frame(lon = c(-65,-65), lat = c(45,47)), size=1, arrow=arrow(length=unit(0.035,'in'), angle=30, type='open'))
low_north <- geom_text(data = data.frame(lon = -65, lat = 48.5, lab = 'N'), aes(label=lab), fontface='bold', size = 3)
scale_bar <- scalebar_latlong(latmin=26, lonmin=-75, h=.3, d=400)
booktheme <- theme(legend.title = element_blank(), legend.position = c(0.25, 0.1), legend.key.width = unit(0.1, 'inches'), legend.key.height = unit(0.05, 'inches'), legend.text = element_text(size = 6))

# Draw maps.
betamap <- draw_fia_map_book(dat = bdmapdat, zvar = 'beta_td_sorensen',  
                        the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
                        colscale = colscalebeta,
                        write_to_file = FALSE, img_h = fig_h - 2, img_w = fig_w/2)

ggsave(file.path(fpfig, 'beta_div_50.png'), betamap + booktheme + ggtitle('Beta-diversity'), height = (fig_h - 2)/3, width = fig_w/2, dpi = 400)

alphamap <- draw_fia_map_book(dat = admapdat, zvar = 'shannon',  
                         the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
                         colscale = colscalealpha,
                         write_to_file = FALSE, img_h = fig_h - 2, img_w = fig_w/2)

ggsave(file.path(fpfig, 'alpha_div_50.png'), alphamap + booktheme + ggtitle('Alpha-diversity'), height = (fig_h - 2)/3, width = fig_w/2, dpi = 400)

gammamap <- draw_fia_map_book(dat = gdmapdat, zvar = 'shannon',  
                         the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
                         colscale = colscalegamma,
                         write_to_file = FALSE, img_h = fig_h - 2, img_w = fig_w/2)

ggsave(file.path(fpfig, 'gamma_div_50.png'), gammamap + booktheme + ggtitle('Gamma-diversity'), height = (fig_h - 2)/3, width = fig_w/2, dpi = 400)

elevmap <- draw_fia_map_book(dat = edmapdat, zvar = 'sd',  
                        the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
                        colscale = colscaleelev,
                        write_to_file = FALSE, img_h = fig_h - 2, img_w = fig_w/2)

ggsave(file.path(fpfig, 'elev_div_50.png'), elevmap + booktheme + ggtitle('Elevation variability'), height = (fig_h - 2)/3, width = fig_w/2, dpi = 400)


# Scatter plots -----------------------------------------------------------

# Plot the relationships.

adplotdat <- ad %>%
  dplyr::filter(east, radius == rad, !is.na(shannon), lon > -100) %>% 
  mutate(shannon = exp(shannon)) %>%
  rename(alpha = shannon) %>%
  dplyr::select(PLT_CN, alpha)
bdplotdat <- bd %>%
  dplyr::filter(east, radius == rad, !is.na(beta_td_sorensen), lon > -100) %>% 
  rename(beta = beta_td_sorensen) %>%
  dplyr::select(PLT_CN, beta)
gdplotdat <- gd %>%
  dplyr::filter(east, radius == rad, !is.na(shannon), lon > -100) %>% 
  mutate(shannon = exp(shannon)) %>%
  rename(gamma = shannon) %>%
  dplyr::select(PLT_CN, gamma)
edplotdat <- ed %>%
  dplyr::filter(east, radius == rad, !is.na(sd), lon > -100, sd <= 500) %>%
  dplyr::select(PLT_CN, sd)
  
scatter_dat <- Reduce(full_join, list(adplotdat, bdplotdat, gdplotdat, edplotdat))

library(reshape2)
scatter_dat <- melt(scatter_dat, id.vars = c('PLT_CN','sd'), variable.name='diversity_type', value.name='diversity_value')

columnscatter <- ggplot(scatter_dat, aes(x = sd, y = diversity_value)) +
  facet_grid(diversity_type ~ ., labeller = labeller(diversity_type = function(x) paste0(x, '-diversity')), scales = 'free_y') +
  geom_hex() +
  scale_fill_gradient(low = 'gray95', high = 'black') +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  scale_x_continuous(name = 'Elevation variability', limits=c(0,400)) +
  scale_y_continuous(name = 'Diversity', expand = c(0,0)) +
  theme(strip.background = element_blank(), legend.position = 'bottom', legend.text = element_text(size=9)) +
  panel_border(colour='black')

ggsave(file.path(fpfig, 'scatter_column.png'), columnscatter, height = fig_h-2, width=fig_w/2, dpi=400)
