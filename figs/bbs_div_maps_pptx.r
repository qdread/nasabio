# Make some bbs diversity maps with the most up to date data and with a cool black background for the pptx
# QDR 01 Feb 2018

# Added 01 Feb 2018: Black map theme.
blackmaptheme_bbs <- theme_bw() + theme(text = element_text(color = 'white'),
                                        panel.grid = element_blank(), 
                                        legend.position = c(0.1, 0.05),
                                        legend.background = element_rect(fill='transparent', color='transparent'),
                                        legend.direction = 'horizontal',
                                        axis.title = element_blank(),
                                        strip.background = element_blank(),
                                        panel.background = element_rect(fill = 'black'),
                                        panel.border = element_rect(color = 'white'),
                                        plot.background = element_rect(fill = 'black'),
                                        plot.title = element_text(color = 'white'),
                                        strip.text = element_text(color = 'white'))

draw_bbs_map <- function(dat, zvar, colscale, by_rad = TRUE,
                         latbds = c(25,50), lonbds = c(-125, -67), 
                         maptheme = whitemaptheme_bbs, 
                         xsc = scale_x_continuous(breaks = seq(-120, -70, by=10), labels = deglab( seq(-120, -70, by=10), 'W')), 
                         ysc = scale_y_continuous(breaks = c(25,35,45), labels = deglab(c(25, 35, 45), 'N')), 
                         maptitle = '',
                         fp = '~', fname = 'plot.png', write_to_file = TRUE, img_h = 7, img_w = 12) {
  the_map <- ggplot(dat, aes(x = lon, y = lat))
  if (by_rad) {
    the_map <- the_map + facet_wrap(~ radius, ncol = 2, labeller = labeller(radius = function(x) paste(as.integer(x)/1000, 'km')))
  }
  the_map <- the_map +
    borders('world', 'canada', fill = 'gray50', color = 'black') +
    borders('world', 'mexico', fill = 'gray50', color = 'black') +
    borders('world', 'usa', fill = 'gray50', color = 'black') +
    borders('state', color = 'black') +
    geom_point(aes_string(color = zvar), size = 0.75) +
    scale_fill_manual(values = c('black', 'white')) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    maptheme + colscale + panel_border(colour = 'black') +
    xsc + ysc + 
    ggtitle(maptitle) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
  if (write_to_file) ggsave(file.path(fp, fname), the_map, height = img_h, width = img_w, dpi = 400) 
  if (!write_to_file) return(the_map)
}


ad %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'richness',  
               colscale = scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, taxonomic, incidence-based (richness)',
               fp = fpfig, fname = 'black_alpha_div_tax.png', img_h = img_height, maptheme = blackmaptheme_bbs)

# By radius
gd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'richness',  
               colscale = scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, taxonomic, incidence-based (richness)',
               fp = fpfig, fname = 'black_gamma_div_tax.png', img_h = img_height, maptheme = blackmaptheme_bbs)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'elevation') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'Elevation\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: elevation',
               fp = fpfig, fname = 'black_geodiv_elev.png', img_h = img_height, maptheme = blackmaptheme_bbs)


bd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'beta_td_pairwise_pa',  
               colscale = scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS beta-diversity, taxonomic, incidence-based',
               fp = fpfig, fname = 'black_beta_div_tax.png', img_h = img_height, maptheme = blackmaptheme_bbs)
