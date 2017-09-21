# Functions for BBS map drawing
# QDR 14 Sep 2017
# Last modified 21 Sep 2017 (ensure 2 columns)

# This will plot maps with the given radius.
# A single row won't work well for BBS, because it's the whole US.
# Instead, we'll make a wrap

deglab <- function(n, d) {
  paste0(as.character(n), '° ', d)
}

draw_bbs_map <- function(dat, zvar, colscale, by_rad = TRUE,
                             scalebar = scalebar_latlong(latmin=24, lonmin=-75, h=.2, d=400), 
                             latbds = c(25,50), lonbds = c(-125, -67), 
                             maptheme = whitemaptheme_bbs, 
                             xsc = scale_x_continuous(breaks = seq(-120, -70, by=10), labels = deglab( seq(-120, -70, by=10), 'W')), 
                             ysc = scale_y_continuous(breaks = c(25,35,45), labels = deglab(c(25, 35, 45), 'N')), 
                             the_arrow = geom_line(data = data.frame(lon = c(-65,-65), lat = c(46,48)), size=1.5, arrow=arrow(length=unit(0.07,'in'), angle=30, type='open')),
                             northlabel = geom_text(data = data.frame(lon = -65, lat = 48.8, lab = 'N'), aes(label=lab), fontface='bold'), 
                             maptitle = '',
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
    geom_text(data=scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 2) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    maptheme + colscale + panel_border(colour = 'black') +
    xsc + ysc + 
    the_arrow + 
    northlabel + 
    ggtitle(maptitle) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
  if (write_to_file) ggsave(file.path(fp, fname), the_map, height = img_h, width = img_w, dpi = 400) 
  if (!write_to_file) return(the_map)
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
    ylab = latmin+3*h,
    text = c(0, d/2, paste(d,'km'))  )
  
  
  list(bar, labs)
}

whitemaptheme_bbs <- theme_bw() + theme(panel.grid = element_blank(), 
                                    legend.position = c(0.1, 0.05),
                                    legend.background = element_rect(fill='transparent', color='transparent'),
                                    legend.direction = 'horizontal',
                                    axis.title = element_blank(),
                                    strip.background = element_blank())