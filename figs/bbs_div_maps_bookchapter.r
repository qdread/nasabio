# BBS maps formatted for the paper.
# Figures taking up an entire page of 8-1/2" by 11" paper with 2.5 cm margins on all sides

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/bookchapter/' # Directory to save images


fig_h <- 11 - 2*(2.5/2.54)
fig_w <- 8.5 - 2*(2.5/2.54)
# Roughly 6.5" wide by 9" high

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
                         fp = '~', fname = 'plot.png', write_to_file = TRUE, img_h = 7, img_w = 12) {
  the_map <- ggplot(dat, aes(x = lon, y = lat))
  if (by_rad) {
    the_map <- the_map + facet_wrap(~ radius, ncol = 2, labeller = labeller(radius = function(x) paste(as.integer(x)/1000, 'km')))
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

draw_bbs_map(dat = beta_td, zvar = 'total',  
             the_arrow = small_low_arrow, northlabel = low_north, scalebar = scale_bar,
             colscale = scale_colour_gradientn(name = 'Taxonomic beta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
             fp = fpfig, fname = 'beta_div_tax_total.png', img_h = fig_h - 2, img_w = fig_w)
