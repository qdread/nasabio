#### New conceptual paper figures (clean version)
# 07 Dec 2017 QDR

# To plot:
# Maps of alpha, beta, gamma, and elevation diversity, one for each of the 5 radii
# Regressions of alpha, beta and gamma versus elevation diversity, one for each of the 5 radii (show confidence bands)
# Plots showing how the R^2 of the regressions change with radius

# Load data (change file path if loading from HPCC)

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_elev_stats_unfuzzed.csv'))
ad <- read.csv(file.path(fp, 'fia_alpha.csv'))
bd <- read.csv(file.path(fp, 'fia_beta.csv'))
gd <- read.csv(file.path(fp, 'fia_gammadiv.csv'))

library(dplyr)
library(cowplot)

load(file.path(fp, 'fiafitplotdat.R'))

# Combine into a single data frame.
biogeo <- ed %>%
  dplyr::select(PLT_CN, STATECD, COUNTYCD, PLOT, radius, sd) %>%
  filter(radius %in% radii) %>%
  rename(elevation_sd = sd) %>%
  left_join(ad %>% 
              dplyr::select(PLT_CN, STATECD, COUNTYCD, PLOT, radius, richness, shannon) %>% 
              rename(alpha_richness = richness, alpha_diversity = shannon) %>%
              filter(radius %in% radii)) %>%
  left_join(bd %>% 
              dplyr::select(PLT_CN, STATECD, COUNTYCD, PLOT, radius, beta_td_pairwise_pa, beta_td_pairwise) %>% 
              rename(beta_richness = beta_td_pairwise_pa, beta_diversity = beta_td_pairwise) %>%
              filter(radius %in% radii)) %>%
  left_join(gd %>% 
              dplyr::select(PLT_CN, STATECD, COUNTYCD, PLOT, radius, richness, shannon) %>% 
              rename(gamma_richness = richness, gamma_diversity = shannon) %>%
              filter(radius %in% radii))


# Customized function to add scale bars and north arrows to maps.
# See http://stackoverflow.com/questions/39067838/parsimonious-way-to-add-north-arrow-and-scale-bar-to-ggmap

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

whitemaptheme <- theme_bw() + theme(panel.grid = element_blank(), 
                                    legend.position = c(0.3, 0.03),
                                    legend.background = element_rect(fill='transparent', color='transparent'),
                                    legend.direction = 'horizontal',
                                    axis.title = element_blank())

fia_xs <- scale_x_continuous(breaks = c(-125, -120, -115), labels = c('125° W', '120° W', '115° W'))
fia_ys <- scale_y_continuous(breaks = c(35,40,45,50), labels = c('35° N', '40° N', '45° N', '50° N'))

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/conceptual_paper'

latbds = c(33,50)
lonbds <- c(-125, -114)

radii <- c(5, 10, 20, 50, 100)


# Add "fuzzed" coordinates.
fiacoords_fuzzed <- read.csv(file.path(fp, 'fia_pnw_coords_fuzzed.csv')) %>% rename(lon=lonfuzz, lat=latfuzz)

bdmapdat_facet <- bd %>% 
  filter(!is.na(beta_td_pairwise_pa), radius %in% radii) %>% 
  left_join(fiacoords_fuzzed) %>%
  arrange(radius, beta_td_pairwise_pa)
admapdat_facet <- ad %>% 
  filter(radius %in% radii, !is.na(shannon)) %>% 
  mutate(shannon = exp(shannon)) %>%
  left_join(fiacoords_fuzzed) %>%
  arrange(radius, shannon)
edmapdat_facet <- ed %>% 
  filter(radius %in% radii, !is.na(sd)) %>% 
  left_join(fiacoords_fuzzed) %>%
  arrange(radius, sd)
gdmapdat_facet <- gd %>% 
  filter(radius %in% radii, !is.na(shannon)) %>% 
  left_join(fiacoords_fuzzed) %>%
  mutate(shannon = exp(shannon)) %>%
  arrange(radius, shannon)

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.2)(9), breaks = c(0,.5,1), limits=c(0,1))
colscalealpha <- scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9))
colscalegamma <- scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = RColorBrewer::brewer.pal(9, 'YlOrRd'))

ptsize <- 0.05 # Points must be very small so that they don't overlap too much.

# Width of figure must be 169 mm
fwidth <- 169

fiamap_bd_facet <- ggplot(bdmapdat_facet, 
                          aes(x = lon, y = lat)) +
  facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  borders('world', 'canada', fill = 'gray90') +
  borders('world', 'usa', fill = 'gray90') +
  borders('state') +
  geom_point(aes(color = beta_td_pairwise_pa), size = ptsize) +
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
  geom_point(aes(color = sd), size = ptsize) +
  geom_rect(data=cbind(westcoast_scalebar[[1]], radius = 100), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 2.5) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscaleelev + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = 'bottom')

fiamap_ad_facet <- ggplot(admapdat_facet, 
                          aes(x = lon, y = lat)) +
  facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  borders('world', 'canada', fill = 'gray90') +
  borders('world', 'usa', fill = 'gray90') +
  borders('state') +
  geom_point(aes(color = shannon), size = ptsize) +
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
  geom_point(aes(color = shannon), size = ptsize) +
  geom_rect(data=cbind(westcoast_scalebar[[1]], radius = 100), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 2.5) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + colscalegamma + panel_border(colour = 'black') +
  fia_xs + fia_ys +
  geom_line(data = data.frame(radius = 100, lon = c(-114.5,-114.5), lat = c(47.6,48.6)), size=0.5, arrow=arrow(length=unit(1, 'mm'), angle=30, type='open')) +
  geom_text(data = data.frame(radius = 100, lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  theme(strip.background = element_blank(), legend.position = 'bottom')

#Extract Legend : code found on stack overflow
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

leftlegtheme <- theme(legend.position = 'left', legend.key.size = unit(5, 'mm'), legend.direction = 'vertical',
                      legend.title = element_text(size = unit(9,'mm')), legend.text = element_text(size = unit(7, 'mm')),
                      axis.text = element_blank(), axis.ticks = element_blank())

nolegtheme <- theme(legend.position = 'none',
                    axis.text = element_blank(), axis.ticks = element_blank())

legend1 <- g_legend(fiamap_bd_facet + leftlegtheme)
legend2 <- g_legend(fiamap_ed_facet + leftlegtheme)
legend3 <- g_legend(fiamap_gd_facet + leftlegtheme + scale_colour_gradientn(name = 'Gamma-\ndiversity', breaks = c(1, 5, 10, 14.9), colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9)))

library(grid)

png(file.path(fpfig, 'legend1.png'), height=fwidth*0.25, width=20, units = 'mm', res=600)
grid.draw(legend1)
dev.off()

png(file.path(fpfig, 'legend2.png'), height=fwidth*0.25, width=20, units = 'mm', res=600)
grid.draw(legend2)
dev.off()

png(file.path(fpfig, 'legend3.png'), height=fwidth*0.25, width=20, units = 'mm', res=600)
grid.draw(legend3)
dev.off()

# Save plots
ggsave(file.path(fpfig, 'SuppFig_alphadiv_maps.png'), fiamap_ad_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'SuppFig_gammadiv_maps.png'), fiamap_gd_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)

ggsave(file.path(fpfig, 'fig2_toprow3_betamap.png'), fiamap_bd_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_middlerow3_elevmap.png'), fiamap_ed_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_toprow_gammamap.png'), fiamap_gd_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)


# Plot the r2 values.
r2_lm_quant$radius_plot <- r2_lm_quant$radius + rep(c(-1,0,1), each = 5)

r2_plot <- ggplot(r2_lm_quant, aes(x = radius_plot, group = interaction(radius, diversity_type), color = diversity_type)) +
  geom_line(aes(y = r2, group = diversity_type), size = 0.25) +
  geom_segment(aes(xend = radius_plot, y = r2_q25, yend = r2_q75), size = 0.5) +
  geom_point(aes(y = r2), size = 1) +
  scale_color_discrete(name = 'Diversity', labels = c('alpha','beta','gamma')) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  theme(axis.text = element_text(size=7.5), axis.title = element_text(size=9))

ggsave(file.path(fpfig, 'fig2_verybottomrow_r2s.png'), r2_plot, width = fwidth * 0.6, height = fwidth * 0.25, units = 'mm', dpi = 600)

# Get rid of values out of range.
xranges <- biogeo %>% group_by(radius) %>% summarize_at(vars(elevation_sd), funs(sdmin = min, sdmax = max), na.rm = TRUE)
pred_val_lm_quant <- pred_val_lm_quant %>%
  ungroup %>%
  left_join(xranges) %>%
  filter(x <= sdmax)

# Median (pseudo) r2s for plotting
r2_lm_quant <- r2_lm_quant %>%
  ungroup %>%
  group_by(diversity_type, radius) %>%
  mutate(r2expr = as.character(eval(substitute(expression(R^2 == x), list(x = round(r2, 2))))))

hexfill <- scale_fill_gradient(low = 'gray90', high = 'black')
hextheme <- theme(strip.background = element_blank())
border <- panel_border(colour = 'black')
radlabel <- labeller(radius = function(x) paste(as.integer(x), 'km'))

alphaplot <- ggplot(biogeo %>% mutate(alpha_diversity=exp(alpha_diversity))) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = alpha_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'alpha_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 1.5) +
  hexfill + hextheme + border + labs(x = 'Elevation variability (m)', y = 'Alpha-diversity')

betaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = beta_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'beta_diversity'),
            aes(x = Inf, y = -Inf, label = r2expr),
            parse = TRUE, hjust = 1.1, vjust = -0.7) +
  scale_x_continuous(labels=seq(0,1200,300), breaks=seq(0,1200,300)) +
  hexfill + hextheme + border + labs(x = 'Elevation variability (m)', y = 'Beta-diversity')

gammaplot <- ggplot(biogeo %>% mutate(gamma_diversity = exp(gamma_diversity))) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = gamma_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'gamma_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 1.5) +
  scale_x_continuous(breaks = c(0,500,1000)) +
  hexfill + hextheme + border + labs(x = 'Elevation variability (m)', y = 'Gamma-diversity')

# Set panel size of bottom row.
source('~/R/setpanelsize.r') # See https://stackoverflow.com/questions/30571198/how-achieve-identical-facet-sizes-and-scales-in-several-multi-facet-ggplot2-grap/30571289#30571289

pset <- set_panel_size(p = betaplot + theme(legend.position = 'none', axis.text = element_text(size = 7.5), axis.title = element_text(size=10.5), strip.background = element_blank(), strip.text.x = element_blank()), 
                       margin = unit(2,'mm'),
                       width = unit(26.5,'mm'),
                       height = unit(fwidth * 0.22,'mm'))
png(file.path(fpfig, 'fig2_bottomrow_fixedsize2.png'), width=fwidth, height=fwidth*0.33, units = 'mm', res = 600)
grid.draw(pset)
dev.off()

pseta <- set_panel_size(p =  alphaplot + theme(legend.position = 'none', axis.text = element_text(size = 7.5), axis.title = element_text(size=10.5), strip.background = element_blank(), strip.text.x = element_text(size=10)), 
                        margin = unit(2,'mm'),
                        width = unit(26.5,'mm'),
                        height = unit(fwidth * 0.22,'mm'))
png(file.path(fpfig, 'SuppFig_alphadiv_fits.png'), width=fwidth, height=fwidth*0.33, units = 'mm', res = 600)
grid.draw(pseta)
dev.off()

psetg <- set_panel_size(p =  gammaplot + theme(legend.position = 'none', axis.text = element_text(size = 7.5), axis.title = element_text(size=10.5), strip.background = element_blank(), strip.text.x = element_text(size=10)), 
                        margin = unit(2,'mm'),
                        width = unit(26.5,'mm'),
                        height = unit(fwidth * 0.23,'mm'))
png(file.path(fpfig, 'SuppFig_gammadiv_fits.png'), width=fwidth, height=fwidth*0.34, units = 'mm', res = 600)
grid.draw(psetg)
dev.off()

# Alternate gamma plot for putting in figure (no column label)
psetg <- set_panel_size(p =  gammaplot + theme(legend.position = 'none', axis.text = element_text(size = 7.5), axis.title = element_text(size=10.5), strip.background = element_blank(), strip.text.x = element_blank()), 
                        margin = unit(2,'mm'),
                        width = unit(26.5,'mm'),
                        height = unit(fwidth * 0.23,'mm'))
png(file.path(fpfig, 'fig2_bottomrow_alternate_gammafits.png'), width=fwidth, height=fwidth*0.34, units = 'mm', res = 600)
grid.draw(psetg)
dev.off()
