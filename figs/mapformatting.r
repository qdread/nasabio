# FIA maps formatted for paper.

# Customized function to add scale bars and north arrows to maps.
# See http://stackoverflow.com/questions/39067838/parsimonious-way-to-add-north-arrow-and-scale-bar-to-ggmap

# great circle distance
# lat1 <- 34
# lat2 <- 34
# lon1 <- -116
# lon2 <- -116 + 1.0848
# 
# lat1 <- lat1*pi/180
# lat2 <- lat2*pi/180
# lon1 <- lon1*pi/180
# lon2 <- lon2*pi/180
# a <- sin((lat1 - lat2)/2)^2 + cos(lat1) * cos(lat2) * sin((lon1 - lon2)/2)^2
# d <- 6371 * 2 * atan2(sqrt(a), sqrt(1-a))

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


westcoast_scalebar <- scalebar_latlong(latmin=33.25, lonmin=-117, h=.2, d=200)

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.3)(9))
colscalealpha <- scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.3)(9))
colscalegamma <- scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.3)(9))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = RColorBrewer::brewer.pal(9, 'YlOrRd'))


whitemaptheme <- theme_bw() + theme(panel.grid = element_blank(), 
                                      legend.position = c(0.3, 0.03),
                                      legend.background = element_rect(fill='transparent', color='transparent'),
                                      legend.direction = 'horizontal',
                                      axis.title = element_blank())

fia_xs <- scale_x_continuous(breaks = c(-125, -120, -115), labels = c('125° W', '120° W', '115° W'))
fia_ys <- scale_y_continuous(breaks = c(35,40,45,50), labels = c('35° N', '40° N', '45° N', '50° N'))

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/betadiv/data_paper_maps'

latbds = c(33,50)
lonbds <- c(-125, -114)

radii <- c(5, 10, 20, 50, 100)

for (rad in radii) {
  
  bdmapdat <- bd %>% filter(radius == rad, !is.na(beta_pairwise_abundance)) %>% arrange(beta_pairwise_abundance)
  admapdat <- ad %>% filter(radius == rad, !is.na(shannon)) %>% arrange(shannon)
  edmapdat <- bd %>% filter(radius == rad, !is.na(sd_elev)) %>% arrange(sd_elev)
  gdmapdat <- gd %>% filter(radius == rad, !is.na(shannon)) %>% arrange(shannon)
  
  fiamap_bd <- ggplot(bdmapdat, 
                      aes(x = lon, y = lat)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(aes(color = beta_pairwise_abundance), size = 0.75) +
    geom_rect(data=westcoast_scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
              show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
    geom_text(data=westcoast_scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    whitemaptheme + colscalebeta + panel_border(colour = 'black') +
    fia_xs + fia_ys +
    geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')) +
    geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('fia_map_',rad,'km_beta.png')
 # ggsave(file.path(fpfig, fname), fiamap_bd, height = 8, width = 6, dpi = 400)
  
  fiamap_ad <- ggplot(admapdat, 
                      aes(x = lon, y = lat)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(aes(color = shannon), size = 0.75) +
    geom_rect(data=westcoast_scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
              show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
    geom_text(data=westcoast_scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    whitemaptheme + colscalealpha + panel_border(colour = 'black') +
    fia_xs + fia_ys +
    geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')) +
    geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('fia_map_',rad,'km_alpha.png')
 # ggsave(file.path(fpfig, fname), fiamap_ad, height = 8, width = 6, dpi = 400)
  
  fiamap_gd <- ggplot(gdmapdat, 
                      aes(x = lon, y = lat)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(aes(color = shannon), size = 0.75) +
    geom_rect(data=westcoast_scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
              show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
    geom_text(data=westcoast_scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    whitemaptheme + colscalegamma + panel_border(colour = 'black') +
    fia_xs + fia_ys +
    geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')) +
    geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('fia_map_',rad,'km_gamma.png')
  ggsave(file.path(fpfig, fname), fiamap_gd, height = 8, width = 6, dpi = 400)
  
  fiamap_ed <- ggplot(edmapdat, 
                      aes(x = lon, y = lat)) +
    borders('world', 'canada', fill = 'gray70') +
    borders('world', 'usa', fill = 'gray70') +
    borders('state') +
    geom_point(aes(color = sd_elev), size = 0.75) +
    geom_rect(data=westcoast_scalebar[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
              show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
    geom_text(data=westcoast_scalebar[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
    coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
    whitemaptheme + colscaleelev + panel_border(colour = 'black') +
    fia_xs + fia_ys +
    geom_line(data = data.frame(lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')) +
    geom_text(data = data.frame(lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
    ggtitle(paste(rad, 'km radius'))
  fname <- paste0('fia_map_',rad,'km_elevdiversity.png')
 # ggsave(file.path(fpfig, fname), fiamap_ed, height = 8, width = 6, dpi = 400)
  
}