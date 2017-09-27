# BBS maps formatted for the paper.
# Figures taking up an entire page of 8-1/2" by 11" paper with 2.5 cm margins on all sides

library(dplyr)
library(cowplot)
library(reshape2)
source('figs/bbs_map_drawing_fns.r')


fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/bookchapter/' # Directory to save images
fpdata <- 'C:/Users/Q/Dropbox/projects/nasabiodiv' # Can be changed to /mnt/research/nasabio/data/bbs to get data from HPCC

# Use the single-year diversities.
bd <- read.csv(file.path(fpdata, 'bbs_betapart_1year.csv'), stringsAsFactors = FALSE)
ad <- read.csv(file.path(fpdata, 'bbs_alpha_1year.csv'), stringsAsFactors = FALSE)
gd <- read.csv(file.path(fpdata, 'bbs_gamma_1year.csv'), stringsAsFactors = FALSE)
ed <- read.csv(file.path(fpdata, 'bbs_geodiversity_stats.csv'), stringsAsFactors = FALSE) # all geodiversity

fig_h <- 11 - 2*(2.5/2.54)
fig_w <- 8.5 - 2*(2.5/2.54)
# Roughly 6.5" wide by 9" high


radii <- c(50000, 75000, 150000)

beta_td <- bd %>%
  filter(diversity == 'taxonomic', radius %in% radii, !is.na(beta)) %>%
  select(-diversity) %>%
  dcast(rteNo + lon + lat + radius ~ partition) %>%
  mutate(prop_nested = nestedness/total,
         prop_replace = replacement/total)

alpha_td <- ad %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii)

gamma_td <- gd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii)

elev_sd <- ed %>%
  mutate(radius = radius * 1000) %>%
  filter(variable == 'elevation', radius %in% radii)

#### Edited map drawing fn
draw_bbs_map <- function(dat, zvar, colscale, by_rad = TRUE,
                         scalebar = scalebar_latlong(latmin=24, lonmin=-75, h=.2, d=400), 
                         latbds = c(25,50), lonbds = c(-125, -67), 
                         maptheme = whitemaptheme_bbs, 
                         xsc = scale_x_continuous(breaks = seq(-120, -70, by=10), labels = deglab( seq(-120, -70, by=10), 'W')), 
                         ysc = scale_y_continuous(breaks = c(25,35,45), labels = deglab(c(25, 35, 45), 'N')), 
                         the_arrow = geom_line(data = data.frame(lon = c(-65,-65), lat = c(46,48)), size=1.5, arrow=arrow(length=unit(0.07,'in'), angle=30, type='open')),
                         northlabel = geom_text(data = data.frame(lon = -65, lat = 48.8, lab = 'N'), aes(label=lab), fontface='bold'), 
                         maptitle = NULL,
                         fp = '~', fname = 'plot.png', write_to_file = TRUE, img_h = 7, img_w = 12, n_columns = 2) {
  the_map <- ggplot(dat, aes(x = lon, y = lat))
  if (by_rad) {
    the_map <- the_map + facet_wrap(~ radius, ncol = n_columns, labeller = labeller(radius = function(x) paste(as.integer(x)/1000, 'km')))
  }
  the_map <- the_map +
    borders('world', 'canada', fill = 'gray90') +
    borders('world', 'mexico', fill = 'gray90') +
    borders('world', 'usa', fill = 'gray90') +
    borders('state') +
    geom_point(aes_string(color = zvar), size = 0.75) +
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
scale_bar <- scalebar_latlong(latmin=23.6, lonmin=-75, h=.3, d=400)

draw_bbs_map(dat = beta_td %>% arrange(total), zvar = 'total',  
             the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
             colscale = scale_colour_gradientn(name = 'Taxonomic beta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
             fp = fpfig, fname = 'beta_div_tax_1column.png', img_h = fig_h - 2, img_w = fig_w/2, n_columns = 1)

draw_bbs_map(dat = alpha_td %>% arrange(richness), zvar = 'richness',  
             the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
             colscale = scale_colour_gradientn(name = 'Taxonomic alpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
             fp = fpfig, fname = 'alpha_div_tax_1column.png', img_h = fig_h - 2, img_w = fig_w/2, n_columns = 1)

draw_bbs_map(dat = gamma_td %>% arrange(richness), zvar = 'richness',  
             the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
             colscale = scale_colour_gradientn(name = 'Taxonomic gamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
             fp = fpfig, fname = 'gamma_div_tax_1column.png', img_h = fig_h - 2, img_w = fig_w/2, n_columns = 1)

draw_bbs_map(dat = elev_sd %>% arrange(sd), zvar = 'sd',  
             the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
             colscale = scale_colour_gradientn(name = 'Elevation\nstandard deviation', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
             fp = fpfig, fname = 'elev_sd_1column.png', img_h = fig_h - 2, img_w = fig_w/2, n_columns = 1)

# Bivariate scatterplot
# 3 columns: radius 50k, 75k, 150k
# 3 rows: alpha, beta, gamma

bbs_xy_plot <- function(xdat, xvar, ydat, yvar, xstat, xlims, xbreaks, ylims, ybreaks, xname, yname, radii) {
  xdat %>%
    filter(variable == xvar) %>%
    right_join(ydat) %>%
    mutate(radius = radius*1000) %>%
    filter(radius %in% radii) %>%
    ggplot(aes_string(x = xstat, y = yvar)) +
    geom_point(alpha = 0.05, size = 0.25) +
    stat_smooth(color = 'red', se = FALSE, method = 'auto') +
    facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(as.integer(x)/1000, 'km'))) +
    scale_x_continuous(name = xname, limits=xlims, breaks=xbreaks) +
    scale_y_continuous(name = yname, limits=ylims, breaks=ybreaks, expand = c(0,0)) +
    theme(strip.background = element_blank()) +
    panel_border(colour='black')
}

alpha_ed_plot <- bbs_xy_plot(xdat = ed %>% select(-richness,-diversity), ydat = ad, xvar = 'elevation', yvar = 'richness', xstat = 'sd', xlims = c(0, 1200), xbreaks = c(0, 500, 1000), ylims = c(0, 132), ybreaks = c(0, 50, 100), xname = 'Elevation standard deviation', yname = 'Taxonomic alpha-diversity', radii = c(50000, 75000, 150000))
beta_ed_plot <- bbs_xy_plot(xdat = ed %>% select(-richness,-diversity), ydat = beta_td %>% mutate(radius=radius/1000), xvar = 'elevation', yvar = 'total', xstat = 'sd', xlims = c(0, 1200), xbreaks = c(0, 500, 1000), ylims = c(0, 1), ybreaks = c(0, .5, 1), xname = 'Elevation standard deviation', yname = 'Taxonomic beta-diversity', radii = c(50000, 75000, 150000))
gamma_ed_plot <- bbs_xy_plot(xdat = ed %>% select(-richness,-diversity), ydat = gd, xvar = 'elevation', yvar = 'richness', xstat = 'sd', xlims = c(0, 1200), xbreaks = c(0, 500, 1000), ylims = c(0, 245), ybreaks = c(0, 50, 100, 150, 200), xname = 'Elevation standard deviation', yname = 'Taxonomic gamma-diversity', radii = c(50000, 75000, 150000))

ggsave(file.path(fpfig, 'plot_alpha_td_by_elev_sd.png'), alpha_ed_plot, height = fig_h/3 - 0.5, width = fig_w, dpi = 400)
ggsave(file.path(fpfig, 'plot_beta_td_by_elev_sd.png'), beta_ed_plot, height = fig_h/3 - 0.5, width = fig_w, dpi = 400)


# Make 3x3 faceted plot with a single dataset for better use of ggplot enabling larger panels.
plotdat <- ed %>%
  filter(variable == 'elevation', radius %in% c(50,75,150)) %>%
  select(rteNo, radius, sd) %>%
  left_join(ad %>% filter(radius %in% c(50,75,150)) %>% select(rteNo, radius, richness) %>% rename(alpha = richness)) %>%
  left_join(gd %>% filter(radius %in% c(50,75,150)) %>% select(rteNo, radius, richness) %>% rename(gamma = richness)) %>%
  left_join(beta_td %>% mutate(radius=radius/1000) %>% filter(radius %in% c(50,75,150)) %>% select(rteNo, radius, total) %>% rename(beta = total)) 
  
plotdat <- melt(plotdat, id.vars = c('rteNo', 'radius', 'sd'), variable.name = 'diversity')

plot3x3 <- ggplot(plotdat %>% mutate(diversity=factor(diversity,levels=c('alpha','beta','gamma'))), aes(x = sd, y = value)) +
  facet_grid(diversity ~ radius, scales = 'free_y', labeller = labeller(radius = function(x) paste(x, 'km'))) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_x_continuous(name = 'Elevation standard deviation', breaks=c(0,500,1000)) +
  scale_y_continuous(name = 'Taxonomic diversity', expand = c(0,0))

ggsave(file.path(fpfig, 'plot_3x3.png'), plot3x3, height = fig_h - 2, width = fig_w, dpi = 400)  
