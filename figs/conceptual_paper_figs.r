#### New conceptual paper figures (clean version)
# 07 Dec 2017 QDR

# Modified 27 Aug 2018: formatting changes requested by GEB reviewer. (also note new file paths to use old PNW data)
# Modified 29 Nov 2018: use newer FIA data with macroplots and plantations removed (now we must manually subset the PNW out)
# Modified 03 Dec 2018: edit axes and color scales for GEB resubmission.

# To plot:
# Maps of alpha, beta, gamma, and elevation diversity, one for each of the 5 radii
# Regressions of alpha, beta and gamma versus elevation diversity, one for each of the 5 radii (show confidence bands)
# Plots showing how the R^2 of the regressions change with radius

# Load data
fp <- '~/Dropbox/projects/nasabiodiv/fia_unfuzzed'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_usa_elev_only.csv'))
ad <- read.csv(file.path(fp, 'fiausa_natural_alpha.csv'))
bd <- read.csv(file.path(fp, 'fiausa_natural_betatd.csv'))
gd <- read.csv(file.path(fp, 'fiausa_natural_gamma.csv'))

library(dplyr)
library(cowplot)

load(file.path(fp, 'fiafitplotdat_pnw.R'))

# Load state codes
fia_statecodes <- read.csv(file.path(fp, 'fiastatecodes.csv'))

radii <- c(5, 10, 20, 50, 100)

# Filter ed dataframe by state code to get only the PNW states, and only the non plantation plots.
pnw_states <- c(6, 41, 53) # Cal, Ore, and Wash.
ed <- ed %>%
  filter(PLT_CN %in% ad$PLT_CN) %>%
  left_join(fia_statecodes) %>%
  filter(STATECD %in% pnw_states) 

# Combine into a single data frame.
biogeo <- ed %>%
  dplyr::select(PLT_CN, STATECD, COUNTYCD, radius, sd) %>%
  filter(radius %in% radii) %>%
  rename(elevation_sd = sd) %>%
  left_join(ad %>% 
              dplyr::select(PLT_CN, radius, richness, shannon) %>% 
              rename(alpha_richness = richness, alpha_diversity = shannon) %>%
              filter(radius %in% radii)) %>%
  left_join(bd %>% 
              dplyr::select(PLT_CN, radius, beta_td_pairwise_pa, beta_td_pairwise) %>% 
              rename(beta_richness = beta_td_pairwise_pa, beta_diversity = beta_td_pairwise) %>%
              filter(radius %in% radii)) %>%
  left_join(gd %>% 
              dplyr::select(PLT_CN, radius, richness, shannon) %>% 
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

fpfig <- '~/google_drive/NASABiodiversityWG/Figures/conceptual_paper'

latbds = c(33,50)
lonbds <- c(-125, -114)

# Add "fuzzed" coordinates.
fiacoords_fuzzed <- read.csv(file.path(fp, 'fia_fuzzed_coords.csv')) %>% rename(lon=lon_fuzzed, lat=lat_fuzzed)

bdmapdat_facet <- biogeo %>% 
  filter(!is.na(beta_diversity), radius %in% radii) %>% 
  left_join(fiacoords_fuzzed) %>%
  rename(beta_td_pairwise = beta_diversity) %>%
  arrange(radius, beta_td_pairwise)
admapdat_facet <- biogeo %>% 
  filter(radius %in% radii, !is.na(alpha_diversity)) %>% 
  mutate(shannon = exp(alpha_diversity)) %>%
  left_join(fiacoords_fuzzed) %>%
  arrange(radius, shannon)
edmapdat_facet <- biogeo %>% 
  filter(radius %in% radii, !is.na(elevation_sd)) %>% 
  rename(sd = elevation_sd) %>%
  left_join(fiacoords_fuzzed) %>%
  arrange(radius, sd)
gdmapdat_facet <- biogeo %>% 
  filter(radius %in% radii, !is.na(gamma_diversity)) %>% 
  left_join(fiacoords_fuzzed) %>%
  mutate(shannon = exp(gamma_diversity)) %>%
  arrange(radius, shannon)

colscalebeta <- scale_colour_gradientn(name = 'Beta-\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.1)(9), breaks = c(0,.5,1), limits=c(0,1))
colscalebeta_transf <- scale_colour_gradientn(name = 'Beta-\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9), labels = c(0, 0.1, 0.5, 0.9, 1), breaks = qlogis(c(0.00001, 0.1, 0.5, 0.9, 0.99999)), limits=qlogis(c(0.00001, 0.99999)))
colscalealpha <- scale_colour_gradientn(name = 'Alpha-\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 2)(9), breaks = c(1, 2, 4, 6), labels = c(1,2,4,6), limits = c(1, 6))
colscalegamma <- scale_colour_gradientn(name = 'Gamma-\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(1, 6, 11, 16), labels=c(1,6,11,16), limits = c(1, 16))
colscaleelev <- scale_colour_gradientn(name = 'Elevation\nvariability', colours = RColorBrewer::brewer.pal(9, 'YlOrRd'), breaks = c(0, 300, 600, 900, 1200), labels = c(0,300,600,900,1200), limits = c(0, 1210))

ptsize <- 0.05 # Points must be very small so that they don't overlap too much.

# Width of figure must be 169 mm
fwidth <- 169

westcoast_scalebar <- scalebar_latlong(latmin=33.25, lonmin=-121, h=.2, d=500)
westcoast_scalebar[[2]] <- transform(westcoast_scalebar[[2]], ylab = ylab + 0.22)

fiamap_bd_facet <- ggplot(bdmapdat_facet %>% filter(beta_td_pairwise <= 0.99999), 
                          aes(x = lon, y = lat)) +
  facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
  borders('world', 'canada', fill = 'gray90') +
  borders('world', 'usa', fill = 'gray90') +
  borders('state') +
  geom_point(aes(color = qlogis(beta_td_pairwise)), size = ptsize) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  whitemaptheme + 
  colscalebeta_transf + 
  panel_border(colour = 'black') +
  fia_xs + fia_ys +
  geom_rect(data=cbind(westcoast_scalebar[[1]], radius = 100), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F,
            show.legend = F,  color = "black", fill = westcoast_scalebar[[1]]$fill.col) +
  geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 1.5) +
  geom_line(data = data.frame(radius = 100, lon = c(-114.5,-114.5), lat = c(47.6,48.6)), size=0.5, arrow=arrow(length=unit(1, 'mm'), angle=30, type='open')) +
  geom_text(data = data.frame(radius = 100, lon = -114.5, lat = 49.2, lab = 'N'), aes(label=lab), fontface='bold') +
  theme(strip.background = element_blank())


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
  geom_text(data=cbind(westcoast_scalebar[[2]], radius = 100), aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F, size = 1.5) +
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
#legend3 <- g_legend(fiamap_gd_facet + leftlegtheme + scale_colour_gradientn(name = 'Gamma-\ndiversity', breaks = c(10, 20, 30, 40), colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 2)(9)))
legend3 <- g_legend(fiamap_gd_facet + leftlegtheme)

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
ggsave(file.path(fpfig, 'SuppFig_betadiv_maps.png'), fiamap_bd_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'SuppFig_gammadiv_maps.png'), fiamap_gd_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)

ggsave(file.path(fpfig, 'fig2_toprow3_betamap.png'), fiamap_bd_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_middlerow3_elevmap.png'), fiamap_ed_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'fig2_toprow_gammamap.png'), fiamap_gd_facet + nolegtheme, width = fwidth - 20, height = fwidth * 0.4, units = 'mm', dpi = 600)

# Save plots as pdf too
ggsave(file.path(fpfig, 'pdfs/SuppFig_alphadiv_maps.pdf'), fiamap_ad_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'pdfs/SuppFig_betadiv_maps.pdf'), fiamap_bd_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'pdfs/SuppFig_gammadiv_maps.pdf'), fiamap_gd_facet + leftlegtheme, width = fwidth, height = fwidth * 0.4, units = 'mm', dpi = 600)

# Plot the r2 values.
r2_lm_quant$radius_plot <- r2_lm_quant$radius + rep(c(-1,0,1), each = 5)

r2_plot <- ggplot(r2_lm_quant, aes(x = radius_plot, group = interaction(radius, diversity_type), color = diversity_type)) +
  geom_line(aes(y = r2, group = diversity_type), size = 0.25) +
  geom_segment(aes(xend = radius_plot, y = r2_q25, yend = r2_q75), size = 0.5) +
  geom_point(aes(y = r2), size = 1) +
  scale_color_discrete(name = 'Diversity', labels = c('alpha','beta','gamma')) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  scale_x_continuous(breaks = radii) +
  theme(axis.text = element_text(size=7.5), axis.title = element_text(size=9))

ggsave(file.path(fpfig, 'fig2_verybottomrow_r2s.png'), r2_plot, width = fwidth * 0.6, height = fwidth * 0.25, units = 'mm', dpi = 600)
ggsave(file.path(fpfig, 'pdfs/SuppFig_r2s.pdf'), r2_plot, width = fwidth * 0.6, height = fwidth * 0.25, units = 'mm', dpi = 600)

# Plot the slopes too.
# Modified 10 Jan 2019: change the shapes of the symbols so it looks OK in B&W
coef_quant$radius_plot <- coef_quant$radius + rep(c(-1,0,1), each = 5)

coef_plot <- ggplot(coef_quant, aes(x = radius_plot, group = interaction(radius, diversity_type), color = diversity_type)) +
  geom_line(aes(y = coef, group = diversity_type), size = 0.25) +
  geom_segment(aes(xend = radius_plot, y = coef_q25, yend = coef_q75), size = 0.5) +
  geom_point(aes(y = coef, shape = diversity_type), size = 1.25) +
  scale_color_discrete(name = 'Diversity', labels = c('alpha','beta','gamma')) +
  scale_shape_discrete(name = 'Diversity', labels = c('alpha','beta','gamma')) +
  labs(x = 'Radius (km)', y = 'Slope coefficient') +
  scale_x_continuous(breaks = radii) +
  theme(axis.text = element_text(size=7.5), axis.title = element_text(size=9))

ggsave(file.path(fpfig, 'fig2_verybottomrow_coefs.png'), coef_plot, width = fwidth * 0.6, height = fwidth * 0.25, units = 'mm', dpi = 600)

# Added 09 May: Edited version of coefficient plot for the proposal.
modified_plot <- coef_plot + 
	labs(x = 'Grain of analysis (radius in km)', y = 'Magnitude of relationship\n(slope coefficient)') +
	theme(legend.position = 'bottom', legend.text = element_text(size=9), legend.title = element_text(size=9))
	
ggsave('~/google_drive/NASABiodiversityWG/Figures/proposal/modified_coef_plot.png', modified_plot, width = 3, height = 3, dpi = 400)	

# Separate coefficient plots for each variable
coef_plot_1x3 <- ggplot(coef_quant, aes(x = radius)) +
  facet_wrap( ~ diversity_type, scales = 'free', labeller = labeller(diversity_type = c(alpha_diversity = 'alpha', beta_diversity = 'beta', gamma_diversity = 'gamma'))) +
  geom_line(aes(y = coef, group = diversity_type), size = 0.25) +
  geom_segment(aes(xend = radius, y = coef_q25, yend = coef_q75), size = 0.5) +
  geom_point(aes(y = coef), size = 1) +
  scale_x_continuous(breaks=radii, labels=radii) +
  labs(x = 'Radius (km)', y = 'Slope coefficient') +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=9), strip.background = element_blank(),
        strip.text.x=element_text(margin=margin(b=5)))

ggsave(file.path(fpfig, 'coefs_separate.png'), coef_plot_1x3, width = fwidth * 0.8, height = fwidth * 0.25, units = 'mm', dpi = 600)

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

alphaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = alpha_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'alpha_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 1.5) +
  scale_x_continuous(labels=seq(0,1200,300), breaks=seq(0,1200,300)) +
  scale_y_continuous(breaks = c(0, .5, 1, 1.5), labels = round(exp(c(0, .5, 1, 1.5)),1)) +
  coord_cartesian(ylim = c(0, 1.5)) +
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

gammaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = gamma_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'gamma_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.1, vjust = 1.5) +
  scale_x_continuous(breaks = c(0,500,1000)) +
  scale_y_continuous(breaks = log(c(1,3,7,20)), labels = c(1, 3, 7, 20)) +
  coord_cartesian(ylim = c(0, 3)) +
  hexfill + hextheme + border + labs(x = 'Elevation variability (m)', y = 'Gamma-diversity')

# Set panel size of bottom row.
source('~/Documents/R/setpanelsize.r') # See https://stackoverflow.com/questions/30571198/how-achieve-identical-facet-sizes-and-scales-in-several-multi-facet-ggplot2-grap/30571289#30571289

library(grid)

pset <- set_panel_size(p = betaplot + theme(legend.position = 'none', axis.text = element_text(size = 7.5), axis.title = element_text(size=10.5), strip.background = element_blank(), strip.text.x = element_blank()), 
                       margin = unit(2,'mm'),
                       width = unit(26.5,'mm'),
                       height = unit(fwidth * 0.22,'mm'))
png(file.path(fpfig, 'SuppFig_betadiv_fits.png'), width=fwidth, height=fwidth*0.33, units = 'mm', res = 600)
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

# Write fit plots as pdfs too
pdf(file.path(fpfig, 'pdfs/SuppFig_betadiv_fits.pdf'), width=fwidth/25.4, height=fwidth*0.33/25.4)
grid.draw(pset)
dev.off()
pdf(file.path(fpfig, 'pdfs/SuppFig_alphadiv_fits.pdf'), width=fwidth/25.4, height=fwidth*0.33/25.4)
grid.draw(pseta)
dev.off()
