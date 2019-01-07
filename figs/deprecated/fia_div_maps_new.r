# Figures with new diversity metrics for FIA plots.
# Made with both abundance-weighted and non-abundance weighted diversity.
# Beta diversity is new, while alpha and gamma are both old but now both abundance-weighted and presence-absence are shown.

# Use the newer background types for consistency.


# Load data ---------------------------------------------------------------

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps/'

bd <- read.csv(file.path(fp, 'fia_betapart_to100.csv'), stringsAsFactors = FALSE) # 17 Aug: only go up to 100km
ed <- read.csv(file.path(fp, 'fia_elev_stats_noalaska.csv'), stringsAsFactors = FALSE)
ad <- read.csv(file.path(fp, 'fia_alpha.csv'), stringsAsFactors = FALSE)
gd <- read.csv(file.path(fp, 'fia_gammadiv.csv'), stringsAsFactors = FALSE)

library(dplyr)

bd <- bd %>%
  mutate(radius = radius/1000) %>% # put radius in km
  filter(radius %in% c(5, 10, 20, 50, 100)) %>% # reduce to more manageable size
  left_join(ed)

ad <- ad %>% left_join(ed)
gd <- gd %>% left_join(ed)


# Methods comparison ------------------------------------------------------

# Before going ahead with maps, we need to compare the methods.
# Compare Sorensen and Jaccard total beta-diversity.
# Also, compare Podani and Baselga partitions, with both Sorensen and Jaccard.

library(cowplot)
library(reshape2)

fsc <- scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'), bias=2)(10))

# Comparison of Sorensen and Jaccard

SJdat <- bd %>% 
  filter(divtype == 'total', family == 'podani') %>%
  select(PLT_CN, radius, index, abundance, beta) %>%
  mutate(abund_name = c('presence_absence', 'abundance_weighted')[abundance+1])

SJdat <- dcast(SJdat, PLT_CN + radius + abund_name ~ index, value.var = 'beta')

ggplot(SJdat, aes(x = jaccard, y = sorensen)) +
  geom_hex() + fsc +
  stat_smooth(method = 'lm') +
  geom_abline(slope=1, intercept=0, color='red', linetype='dotted') +
  facet_grid(abund_name ~ radius) +
  theme_bw() +
  ggtitle('Comparison of Sorensen and Jaccard indices')

# Comparison of Baselga and Podani decompositions

BPdat <- bd %>%
  filter(divtype %in% c('replacement_proportion')) %>%
  select(PLT_CN, radius, family, index, divtype, abundance, beta) %>%
  mutate(abund_name = c('presence_absence', 'abundance_weighted')[abundance+1])

BPdat_sorensen <- filter(BPdat, index == 'sorensen')
BPdat_jaccard <- filter(BPdat, index == 'jaccard')

BPdat_sorensen <- dcast(BPdat_sorensen, PLT_CN + radius + abund_name ~ family, value.var = 'beta')
BPdat_jaccard <- dcast(BPdat_jaccard, PLT_CN + radius + abund_name ~ family, value.var = 'beta')

ggplot(BPdat_sorensen, aes(x = baselga, y = podani)) +
  geom_hex() + fsc +
  stat_smooth(method = 'lm') +
  geom_abline(slope=1, intercept=0, color='red', linetype='dotted') +
  facet_grid(abund_name ~ radius) +
  theme_bw() +
  ggtitle('Comparison of Baselga and Podani partitioning: proportion due to species replacement', subtitle = 'Sorensen index')

ggplot(BPdat_jaccard, aes(x = baselga, y = podani)) +
  geom_hex() + fsc +
  stat_smooth(method = 'lm') +
  geom_abline(slope=1, intercept=0, color='red', linetype='dotted') +
  facet_grid(abund_name ~ radius) +
  theme_bw() +
  ggtitle('Comparison of Baselga and Podani partitioning: proportion due to species replacement', subtitle = 'Jaccard index')

# Comparison of abundance-weighted and presence-absence diversity

abunddat <- bd %>%
  filter(family == 'podani', divtype == 'total') %>%
  select(PLT_CN, radius, index, abundance, beta) %>%
  mutate(abund_name = c('presence_absence', 'abundance_weighted')[abundance+1])

abunddat <- dcast(abunddat, PLT_CN + radius + index ~ abund_name, value.var = 'beta')

ggplot(abunddat, aes(x = presence_absence, y = abundance_weighted)) +
  geom_hex() + fsc +
  stat_smooth(method = 'lm') +
  geom_abline(slope=1, intercept=0, color='red', linetype='dotted') +
  facet_grid(index ~ radius) +
  theme_bw() +
  ggtitle('Comparison of presence-based and abundance-based beta-diversity')

# Plotting functions ------------------------------------------------------

# This will plot just some individual maps, rather than the big row of maps formatted for the paper.
# Refer back to mapformatting.r for drawing faceted maps.

draw_fia_map <- function(dat, zvar, rad, colscale, 
                         scalebar = scalebar_latlong(latmin=33.25, lonmin=-117, h=.2, d=200), 
                         latbds = c(33,50), lonbds = c(-125, -114), 
                         maptheme = whitemaptheme, 
                         xsc = scale_x_continuous(breaks = c(-125, -120, -115), labels = c('125° W', '120° W', '115° W')), 
                         ysc = scale_y_continuous(breaks = c(35,40,45,50), labels = c('35° N', '40° N', '45° N', '50° N')), 
                         the_arrow = geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')),
                         northlabel = geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold'), 
                         maptitle = ggtitle(paste(rad, 'km radius')),
                         fp, fname) {
  the_map <- ggplot(dat, 
                      aes(x = lon, y = lat)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(aes_string(color = zvar), size = 0.75) +
    geom_rect(data=scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
              show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
    geom_text(data=scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    maptheme + colscale + panel_border(colour = 'black') +
    xsc + ysc + the_arrow + northlabel + maptitle
  ggsave(file.path(fp, fname), the_map, height = 8, width = 6, dpi = 400)
}

# This will plot a row of maps with the given radius.
draw_fia_map_row <- function(dat, zvar, rad, colscale, 
                             scalebar = scalebar_latlong(latmin=33.25, lonmin=-117, h=.2, d=200), 
                             latbds = c(33,50), lonbds = c(-125, -114), 
                             maptheme = whitemaptheme, 
                             xsc = scale_x_continuous(breaks = c(-125, -120, -115), labels = c('125° W', '120° W', '115° W')), 
                             ysc = scale_y_continuous(breaks = c(35,40,45,50), labels = c('35° N', '40° N', '45° N', '50° N')), 
                             the_arrow = geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')),
                             northlabel = geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold'), 
                             maptitle = '',
                             fp, fname) {
  the_map <- ggplot(dat, aes(x = lon, y = lat)) +
    facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
    borders('world', 'canada', fill = 'gray90') +
    borders('world', 'usa', fill = 'gray90') +
    borders('state') +
    geom_point(aes_string(color = zvar), size = 0.75) +
    geom_rect(data=scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
              show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
    geom_text(data=scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    maptheme + colscale + panel_border(colour = 'black') +
    xsc + ysc + 
    the_arrow + 
    northlabel + 
    ggtitle(maptitle)
  ggsave(file.path(fp, fname), the_map, height = 8, width = 20, dpi = 400) 
}

# function to find delta longitude for two points separated by x km at a certain latitude
deltalong <- function(d, lat) {
  R <- 6371 # earth's radius
  r <- R * cos(lat * pi/180)
  circ <- 2*pi*r # this makes 360 degrees of longitude at that point
  d * 360/circ
}

scalebar_latlong = function(latmin, lonmin, d, h){
  lonmax <- lonmin + deltalong(d, latmin)
  # lonmin,latmin = lower left coordinate of bar
  # h = height of bar
  # d = distance of bar
  
  bar <- data.frame(
    xmin = c(lonmin, (lonmin+lonmax)/2),
    xmax = c((lonmin+lonmax)/2, lonmax),
    ymin = c(latmin, latmin),
    ymax = c(latmin+h, latmin+h),
    z = c(1,0),
    fill.col = c('black','white')
  )
  labs <- data.frame(
    xlab = c(lonmin, (lonmin+lonmax)/2, lonmax),
    ylab = latmin+2*h,
    text = c(0, d/2, paste(d,'km'))  )
  
  
  list(bar, labs)
}

# Plot settings -----------------------------------------------------------

whitemaptheme <- theme_bw() + theme(panel.grid = element_blank(), 
                                    legend.position = c(0.3, 0.03),
                                    legend.background = element_rect(fill='transparent', color='transparent'),
                                    legend.direction = 'horizontal',
                                    axis.title = element_blank())

fia_xs <- scale_x_continuous(breaks = c(-125, -120, -115), labels = c('125° W', '120° W', '115° W'))
fia_ys <- scale_y_continuous(breaks = c(35,40,45,50), labels = c('35° N', '40° N', '45° N', '50° N'))

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.2)(9), breaks = c(0,.5,1), limits=c(0,1))
colscalealpha <- scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9))
colscalegamma <- scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = RColorBrewer::brewer.pal(9, 'YlOrRd'))

westcoast_scalebar <- scalebar_latlong(latmin=33.25, lonmin=-121, h=.2, d=500)
westcoast_scalebar[[2]] <- transform(westcoast_scalebar[[2]], ylab = ylab + 0.22)

# Draw plots --------------------------------------------------------------

radii <- c(5, 10, 20, 50, 100)

bd_incidence <- bd %>% filter(family == 'baselga', index == 'sorensen', !abundance, !is.na(beta), radius %in% radii)
bd_incidence_total <- bd_incidence %>% filter(divtype == 'total') %>% arrange(radius, beta)
bd_incidence_replacement <- bd_incidence %>% filter(divtype == 'replacement_proportion') %>% arrange(radius, beta)
bd_incidence_nested <- bd_incidence %>% filter(divtype == 'nestedness_proportion') %>% arrange(radius, beta)

bd_abundance <- bd %>% filter(family == 'baselga', index == 'sorensen', !abundance, !is.na(beta), radius %in% radii)
bd_abundance_total <- bd_abundance %>% filter(divtype == 'total') %>% arrange(radius, beta)
bd_abundance_replacement <- bd_abundance %>% filter(divtype == 'replacement_proportion') %>% arrange(radius, beta)
bd_abundance_nested <- bd_abundance %>% filter(divtype == 'nestedness_proportion') %>% arrange(radius, beta)

ad_incidence <- ad %>% filter(radius %in% radii, !is.na(richness)) %>% arrange(radius, richness)
ad_abundance <- ad %>% filter(radius %in% radii, !is.na(shannon)) %>% arrange(radius, shannon)

gd_incidence <- gd %>% filter(radius %in% radii, !is.na(richness)) %>% arrange(radius, richness)
gd_abundance <- gd %>% filter(radius %in% radii, !is.na(shannon)) %>% arrange(radius, shannon)

draw_fia_map_row(dat = bd_incidence_total, zvar = 'beta', rad = radii, 
                 colscale = colscalebeta, 
                 maptitle = 'FIA Beta-diversity, incidence-based, total',
                 fp = fpfig, fname = 'beta_div_incidence_total.png')

