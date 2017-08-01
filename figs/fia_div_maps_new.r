# Figures with new diversity metrics for FIA plots.
# Made with both abundance-weighted and non-abundance weighted diversity.
# Beta diversity is new, while alpha and gamma are both old but now both abundance-weighted and presence-absence are shown.

# Use the newer background types for consistency.


# Load data ---------------------------------------------------------------

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_diversity_maps/'

bd <- read.csv(file.path(fp, 'fia_betapart_to150.csv'), stringsAsFactors = FALSE) # Slow.
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

# Comparison of Sorensen and Jaccard

SJdat <- bd %>% 
  filter(divtype == 'total', family == 'podani') %>%
  select(PLT_CN, radius, index, abundance, beta) %>%
  mutate(abund_name = c('presence_absence', 'abundance_weighted')[abundance+1])

SJdat <- dcast(SJdat, PLT_CN + radius + abund_name ~ index, value.var = 'beta')

ggplot(SJdat, aes(x = jaccard, y = sorensen)) +
  geom_point() +
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
  geom_point() +
  stat_smooth(method = 'lm') +
  geom_abline(slope=1, intercept=0, color='red', linetype='dotted') +
  facet_grid(abund_name ~ radius) +
  theme_bw() +
  ggtitle('Comparison of Baselga and Podani partitioning: proportion due to species replacement', subtitle = 'Sorensen index')

ggplot(BPdat_jaccard, aes(x = baselga, y = podani)) +
  geom_point() +
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
  geom_point() +
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
                         lonbds = c(33,50), latbds = c(-125, -114), 
                         maptheme = whitemaptheme, 
                         xsc = scale_x_continuous(breaks = c(-125, -120, -115), labels = c('125° W', '120° W', '115° W')), 
                         ysc = scale_y_continuous(breaks = c(35,40,45,50), labels = c('35° N', '40° N', '45° N', '50° N')), 
                         arrow = geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')),
                         northlabel = geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold'), 
                         maptitle = ggtitle(paste(rad, 'km radius')),
                         fpfig, fname) {
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
    xsc + ysc + arrow + northlabel + maptitle
  ggsave(file.path(fpfig, fname), the_map, height = 8, width = 6, dpi = 400)
}
