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
  #ggsave(file.path(fpfig, fname), fiamap_gd, height = 8, width = 6, dpi = 400)
  
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


########################################################

# Edit 27 June.
# Make the figures according to specific, newer guidelines.
# One big figure with a lot of panels.

bdmapdat_facet <- bd %>% filter(!is.na(beta_pairwise_abundance), radius %in% radii) %>% arrange(radius, beta_pairwise_abundance)
admapdat_facet <- ad %>% filter(radius %in% radii, !is.na(shannon)) %>% arrange(radius, shannon)
edmapdat_facet <- bd %>% filter(radius %in% radii, !is.na(sd_elev)) %>% arrange(radius, sd_elev)
gdmapdat_facet <- gd %>% filter(radius %in% radii, !is.na(shannon)) %>% arrange(radius, shannon)

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.2)(9), breaks = c(0,.5,1), limits=c(0,1))
colscalealpha <- scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9))
colscalegamma <- scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = RColorBrewer::brewer.pal(9, 'YlOrRd'))


fiamap_bd_facet <- ggplot(bdmapdat_facet, 
                    aes(x = lon, y = lat)) +
  facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  borders('world', 'canada', fill = 'gray90') +
  borders('world', 'usa', fill = 'gray90') +
  borders('state') +
  geom_point(aes(color = beta_pairwise_abundance), size = 0.5) +
  #geom_rect(data=cbind(westcoast_scalebar[[1]], radius = 100), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
  #          show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  #geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscalebeta + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  geom_line(data = data.frame(radius = 100, lon = c(-114.5,-114.5), lat = c(47.6,48.6)), size=0.5, arrow=arrow(length=unit(1, 'mm'), angle=30, type='open')) +
  geom_text(data = data.frame(radius = 100, lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  theme(strip.background = element_blank())

westcoast_scalebar <- scalebar_latlong(latmin=33.25, lonmin=-121, h=.2, d=500)
westcoast_scalebar[[2]] <- transform(westcoast_scalebar[[2]], ylab = ylab + 0.22)

fiamap_ed_facet <- ggplot(edmapdat_facet, 
                          aes(x = lon, y = lat)) +
  facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  borders('world', 'canada', fill = 'gray90') +
  borders('world', 'usa', fill = 'gray90') +
  borders('state') +
  geom_point(aes(color = sd_elev), size = 0.5) +
  geom_rect(data=cbind(westcoast_scalebar[[1]], radius = 100), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 2.5) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscaleelev + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  #geom_line(data = data.frame(radius = 100, lon = c(-114.5,-114.5), lat = c(48,49)), size=1.5, arrow=arrow(length=unit(0.1,'in'), angle=30, type='open')) +
  #geom_text(data = data.frame(radius = 100, lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = 'bottom')

fiamap_ad_facet <- ggplot(admapdat_facet, 
                          aes(x = lon, y = lat)) +
  facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  borders('world', 'canada', fill = 'gray90') +
  borders('world', 'usa', fill = 'gray90') +
  borders('state') +
  geom_point(aes(color = shannon), size = 0.5) +
  geom_rect(data=cbind(westcoast_scalebar[[1]], radius = 100), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 2.5) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscalealpha + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  geom_line(data = data.frame(radius = 100, lon = c(-114.5,-114.5), lat = c(47.6,48.6)), size=0.5, arrow=arrow(length=unit(1, 'mm'), angle=30, type='open')) +
  geom_text(data = data.frame(radius = 100, lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  theme(strip.background = element_blank(), legend.position = 'bottom')

fiamap_gd_facet <- ggplot(gdmapdat_facet, 
                          aes(x = lon, y = lat)) +
  facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  borders('world', 'canada', fill = 'gray90') +
  borders('world', 'usa', fill = 'gray90') +
  borders('state') +
  geom_point(aes(color = shannon), size = 0.5) +
  geom_rect(data=cbind(westcoast_scalebar[[1]], radius = 100), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 2.5) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscalegamma + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  geom_line(data = data.frame(radius = 100, lon = c(-114.5,-114.5), lat = c(47.6,48.6)), size=0.5, arrow=arrow(length=unit(1, 'mm'), angle=30, type='open')) +
  geom_text(data = data.frame(radius = 100, lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  theme(strip.background = element_blank(), legend.position = 'bottom')

fia_beta_plot <- ggplot(bd %>% subset(radius %in% radii), aes(x = sd_elev, y = beta_pairwise_abundance)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Taxonomic beta-diversity', expand = c(0,0), limits = c(0,1)) +
  labs(x = 'Elevation variability')

# Width of figure must be 169 mm
fwidth <- 169

fifteenplots <- plot_grid(fiamap_bd_facet, fiamap_ed_facet, fia_beta_plot, nrow = 3, align = 'hv')
ggsave(file.path(fpfig, 'fig2_15plots.png'), fifteenplots, width = fwidth, height = fwidth * 1.2, units = 'mm', dpi = 400)

goodtheme <- theme(legend.position = c(0.153,0.12), legend.key.size = unit(5, 'mm'), legend.box.background = element_rect(fill = 'white'),
                   legend.title = element_text(size = unit(9,'mm')), legend.text = element_text(size = unit(7, 'mm')),
                   axis.text = element_blank(), axis.ticks = element_blank())

leftlegtheme <- theme(legend.position = 'left', legend.key.size = unit(5, 'mm'), legend.direction = 'vertical',
                   legend.title = element_text(size = unit(9,'mm')), legend.text = element_text(size = unit(7, 'mm')),
                   axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(fpfig, 'fig2_toprow_betamap.png'), fiamap_bd_facet + goodtheme, width = fwidth - 10, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_middlerow_elevmap.png'), fiamap_ed_facet + goodtheme, width = fwidth - 10, height = fwidth * 0.4, units = 'mm', dpi = 600)

ggsave(file.path(fpfig, 'fig2_toprow2_betamap.png'), fiamap_bd_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_middlerow2_elevmap.png'), fiamap_ed_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_bottomrow_fits.png'), fia_beta_plot + theme(axis.text = element_text(size = 9), axis.title = element_text(size=10.5)), width = fwidth, height = fwidth * 0.32, units = 'mm', dpi = 600)

# other three panels
pgamfia <- ggplot(fiafitdat, aes(x=radius, y=rsq)) + 
  stat_smooth(method = lm, formula = y ~ x + I(x^2), se=FALSE, size=1, color='red') +
  geom_point(size = 2) + 
  geom_text(data=data.frame(radius=23, rsq=c(.50,.46), diversity=rep(c('alpha','beta','gamma'),each=2), lab=c('alpha-diversity','R^2 == 0.944','beta-diversity','R^2 == .995','gamma-diversity','R^2 == 0.989')), aes(label=lab), parse=TRUE, size = 2.5) +
  facet_wrap(~ diversity) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank(), strip.text = element_blank(), axis.text = element_text(size=7), axis.title = element_text(size=11)) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  scale_x_continuous(breaks=radii) 
ggsave(file.path(fpfig, 'fig3_verybottomrow_gamfits.png'), pgamfia, width = fwidth * 0.75, height = fwidth * 0.25, units = 'mm', dpi = 600)


#####

#edited 28 june. make legend into a separate thing.
#Extract Legend : code found on stack overflow
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend1 <- g_legend(fiamap_bd_facet + leftlegtheme)
legend2 <- g_legend(fiamap_ed_facet + leftlegtheme)

library(grid)

png(file.path(fpfig, 'legend1.png'), height=fwidth*0.25, width=20, units = 'mm', res=600)
grid.draw(legend1)
dev.off()

png(file.path(fpfig, 'legend2.png'), height=fwidth*0.25, width=20, units = 'mm', res=600)
grid.draw(legend2)
dev.off()

# Make plots with no legend.

nolegtheme <- theme(legend.position = 'none',
                      axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(fpfig, 'fig2_toprow3_betamap.png'), fiamap_bd_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_middlerow3_elevmap.png'), fiamap_ed_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_bottomrow_fits.png'), fia_beta_plot + theme(axis.text = element_text(size = 9), axis.title = element_text(size=10.5)), width = fwidth - 10, height = fwidth * 0.32, units = 'mm', dpi = 600)


bplot <- ggplot_build(fia_beta_plot + theme(axis.text = element_text(size = 9), axis.title = element_text(size=10.5)))
btable <- ggplot_gtable(bplot)

source('~/R/setpanelsize.r') # See https://stackoverflow.com/questions/30571198/how-achieve-identical-facet-sizes-and-scales-in-several-multi-facet-ggplot2-grap/30571289#30571289

pset <- set_panel_size(p = fia_beta_plot + theme(axis.text = element_text(size = 9), axis.title = element_text(size=10.5)), 
               margin = unit(2,'mm'),
               width = unit(26.5,'mm'),
               height = unit(fwidth * 0.28,'mm'))
png(file.path(fpfig, 'fig2_bottomrow_fixedsize.png'), width=fwidth, height=fwidth*0.4, units = 'mm', res = 600)
grid.draw(pset)
dev.off()


####################################

# Set panel sizes for each of the three rows and 