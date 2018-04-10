# IALE Figures

# Outline of talk based on figures used

# Methods section
# 1. Schematic illustrating the basic metrics
# 2. Schematic depicting kernels (scale)
# 3. Show how BBS and FIA plots are arranged within the circles, and how we get diversity from them (a, b, g diversity)
# (here some maps of BBS and FIA plots)

# Results section
# 4. Result: comparison of same metric at different scales for TRI and elevational SD
# 5. Result: comparison of TRI versus elevational SD at a number of scales
# 6. Which one explains more variation at the different scales? (has different slope) Is it different for BBS vs FIA and is it different for each diversity type?


# Load data on remote and save locally ------------------------------------

# This section edited on 05 April to use precipitation instead of temperature.
fpfia <- '/mnt/research/nasabio/data/fia/geodiv'

library(dplyr)

fia_elev_long <- read.csv(file.path(fpfia, 'fia_usa_elev.csv'), stringsAsFactors = FALSE) %>%
  select(PLT_CN, variable, radius, mean, sd, min, max) %>%
  filter(grepl('elevation_5k', variable))

fia_other <- read.csv(file.path(fpfia, 'fia_usa_other.csv'), stringsAsFactors = FALSE) %>%
  select(PLT_CN, variable, radius, mean, sd, min, max) %>%
  filter(grepl('human_footprint_5k', variable))

fia_bio5k <- read.csv(file.path(fpfia, 'fia_usa_bio5k.csv'), stringsAsFactors = FALSE) %>%
  select(PLT_CN, variable, radius, mean, sd, min, max) %>%
  filter(grepl('bio12_5k', variable))

fia_vars <- rbind(fia_elev_long, fia_other, fia_bio5k) 
write.csv(fia_vars, file = '/mnt/research/nasabio/temp/fia_vars_iale.csv', row.names = FALSE)

fpbbs <- '/mnt/research/nasabio/data/bbs'

bbs_vars <- read.csv(file.path(fpbbs, 'bbs_geodiversity.csv'), stringsAsFactors = FALSE)

# Keep only elevation, human footprint, and mean annual precipitation.

bbs_vars <- bbs_vars %>%
  filter(grepl('elevation_5k|human_footprint_5k|bio12_5k', variable))
write.csv(bbs_vars, file = '/mnt/research/nasabio/temp/bbs_vars_iale.csv', row.names = FALSE)


# Load data and reshape for plotting if needed ----------------------------

library(dplyr)
library(tidyr)
library(reshape2)

# Use BBS only at first (smaller)
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
bbs_vars <- read.csv(file.path(fp, 'bbs_vars_iale.csv'), stringsAsFactors = FALSE)

bbs_vars <- bbs_vars %>%
  mutate(variable = gsub('human_','human',variable)) %>%
  separate(variable, into = c('variable','resolution','metric')) %>%
  mutate(metric = if_else(is.na(metric), 'raw', metric)) %>%
  mutate_at(c('mean','sd','min','max'), funs(ifelse(variable == 'bio12', ./10, .))) %>%
  select(-richness_geodiv, -diversity_geodiv, -mode, -lon, -lat, -lon_aea, -lat_aea)

bbs_vars_melt <- melt(bbs_vars, id.vars = 1:5, variable.name = 'summary_stat', value.name = 'value')
bbs_vars_cast <- dcast(bbs_vars_melt, rteNo + variable + resolution + radius ~ metric + summary_stat)
bbs_vars_cast_wide <- dcast(bbs_vars_melt, rteNo + variable + resolution ~ metric + summary_stat + radius)

# FIA (a lot bigger, csv is 1gb)
fia_vars <- read.csv(file.path(fp, 'fia_vars_iale.csv'), stringsAsFactors = FALSE)

fia_vars <- fia_vars %>%
  mutate(variable = gsub('human_','human',variable)) %>%
  separate(variable, into = c('variable','resolution','metric')) %>%
  mutate(metric = if_else(is.na(metric), 'raw', metric)) %>%
  mutate_at(c('mean','sd','min','max'), funs(ifelse(variable == 'bio12', ./10, .)))


fia_vars_melt <- melt(fia_vars, id.vars = 1:5, variable.name = 'summary_stat', value.name = 'value')
fia_vars_cast <- dcast(fia_vars_melt, PLT_CN + variable + resolution + radius ~ metric + summary_stat)
fia_vars_cast_wide <- dcast(fia_vars_melt, PLT_CN + variable + resolution ~ metric + summary_stat + radius)

rm(fia_vars, fia_vars_melt)

# Make plots comparing the metrics at different radii ---------------------

library(ggplot2)

radii <- c(5, 10, 20, 50, 100, 200)
geo_vars <- c('bio12', 'elevation', 'humanfootprint')
geo_var_names <- c('mean annual precipitation', 'elevation', 'human footprint index')

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/iale_figs'

threebreaks <- function(limits) {
  x <- max(limits)*0.9
  ndigx <- floor(log10(x))
  lim <- floor(x/(10^ndigx)) * 10^ndigx
    
  c(0, lim/2, lim)
}
scale_x <- scale_x_continuous(breaks = threebreaks)

th_hexplot <- theme_bw() + theme(strip.background = element_rect(fill = NA),
                                 legend.position = 'none')
p_height <- 2.2
p_width <- 6.5

for (i in 1:length(geo_vars)) {

p1 <- fia_vars_cast %>%
  filter(variable == geo_vars[i], radius %in% radii) %>%
  ggplot(aes(x = raw_sd, y = roughness_mean)) +
    facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km radius'))) +
    geom_hex() +
    scale_fill_gradient(low = 'gray80', high = 'black') +
    #theme_bw() +
    #theme(strip.background = element_blank()) +
    #ggtitle('Roughness vs. standard deviation', geo_var_names[i]) +
    labs(x = 'Standard deviation', y = 'Roughness') +
    scale_x + th_hexplot
 
p2 <- fia_vars_cast %>%
  filter(variable == geo_vars[i], radius %in% radii) %>%
  ggplot(aes(x = raw_sd, y = tri_mean)) +
    facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km radius'))) +
    geom_hex() +
    scale_fill_gradient(low = 'gray80', high = 'black') +
    #theme_bw() +
    #theme(strip.background = element_blank()) +
    #ggtitle('Terrain ruggedness index vs. standard deviation', geo_var_names[i]) +
    labs(x = 'Standard deviation', y = 'TRI') +
    scale_x + th_hexplot

p3 <- fia_vars_cast %>%
  filter(variable == geo_vars[i], radius %in% radii) %>%
  ggplot(aes(x = roughness_mean, y = tri_mean)) +
    facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km radius'))) +
    geom_hex() +
    scale_fill_gradient(low = 'gray80', high = 'black') +
    #theme_bw() +
    #theme(strip.background = element_blank()) +
    #ggtitle('Terrain ruggedness index vs. roughness', geo_var_names[i]) +
    labs(x = 'Roughness', y = 'TRI') +
    scale_x + th_hexplot

ggsave(file.path(fpfig, paste('newsize_compare', geo_vars[i], 'roughness_stdev.png', sep = '_')), p1, height = p_height, width = p_width, dpi = 400)
ggsave(file.path(fpfig, paste('newsize_compare', geo_vars[i], 'ruggedness_stdev.png', sep = '_')), p2, height = p_height, width = p_width, dpi = 400)
ggsave(file.path(fpfig, paste('newsize_compare', geo_vars[i], 'ruggedness_roughness.png', sep = '_')), p3, height = p_height, width = p_width, dpi = 400)

}

# Make plot showing how a single metric changes with radius ---------------

### Line plots
set.seed(555)
plot_sample <- sample(unique(fia_vars_cast_wide$PLT_CN), 500)

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii, PLT_CN %in% plot_sample) %>%
  ggplot(aes(x = radius, y = raw_sd, group = PLT_CN)) +
    geom_line(size = 0.1, alpha = 0.2) +
    scale_y_continuous(name = 'Standard deviation', limits = c(-10, 1200), expand = c(0,0)) +
    scale_x_continuous(name = 'Radius (km)', breaks = radii, labels = radii) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggtitle('Standard deviation change with radius', 'Elevation')

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii, PLT_CN %in% plot_sample) %>%
  ggplot(aes(x = radius, y = roughness_mean, group = PLT_CN)) +
    geom_line(size = 0.1, alpha = 0.2) +
    scale_y_continuous(name = 'Roughness', limits = c(-10, 2000), expand = c(0,0)) +
    scale_x_continuous(name = 'Radius (km)', breaks = radii, labels = radii) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggtitle('Roughness change with radius', 'Elevation')    

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii, PLT_CN %in% plot_sample) %>%
  ggplot(aes(x = radius, y = tri_mean, group = PLT_CN)) +
  geom_line(size = 0.1, alpha = 0.2) +
  scale_y_continuous(name = 'TRI', limits = c(-10, 800), expand = c(0,0)) +
  scale_x_continuous(name = 'Radius (km)', breaks = radii, labels = radii) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle('Terrain ruggedness index change with radius', 'Elevation')    

### Individual density plots
for (i in 1:length(geo_vars)) {

p1 <- fia_vars_cast %>%
  filter(variable == geo_vars[i], radius %in% radii) %>%
  ggplot(aes(x = raw_sd, group = radius, fill = factor(radius))) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(name = 'Standard deviation') +
    scale_fill_brewer(name = 'Radius (km)', palette = 'Set1') +
    theme_bw() +
    ggtitle('Standard deviation change with radius', geo_var_names[i])

p2 <- fia_vars_cast %>%
  filter(variable == geo_vars[i], radius %in% radii) %>%
  ggplot(aes(x = roughness_mean, group = radius, fill = factor(radius))) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(name = 'Roughness') +
    scale_fill_brewer(name = 'Radius (km)', palette = 'Set1') +
    theme_bw() +
    ggtitle('Roughness change with radius', geo_var_names[i])

p3 <- fia_vars_cast %>%
  filter(variable == geo_vars[i], radius %in% radii) %>%
  ggplot(aes(x = tri_mean, group = radius, fill = factor(radius))) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(name = 'TRI') +
    scale_fill_brewer(name = 'Radius (km)', palette = 'Set1') +
    theme_bw() +
    ggtitle('Terrain ruggedness index change with radius', geo_var_names[i])

ggsave(file.path(fpfig, paste('density', geo_vars[i], 'stdev.png', sep = '_')), p1, height = 5, width = 5, dpi = 400)
ggsave(file.path(fpfig, paste('density', geo_vars[i], 'roughness.png', sep = '_')), p2, height = 5, width = 5, dpi = 400)
ggsave(file.path(fpfig, paste('density', geo_vars[i], 'ruggedness.png', sep = '_')), p3, height = 5, width = 5, dpi = 400)

}

### Density plot combined in single faceted figure.
# Put mean,sd,min,max into one column
# This doesn't work very well so probably better to combine the plots with cowplot.
metric_names <- c(raw_sd = 'standard deviation', roughness_mean = 'roughness', tri_mean = 'TRI')
variable_names <- c(bio1 = 'temperature', elevation = 'elevation', humanfootprint = 'human footprint', bio12 = 'precipitation')

p_dens_all <- 
  melt(fia_vars, id.vars = 1:5, measure.vars = 6:9, variable.name = 'summary_stat') %>% 
  mutate(metric_stat = paste(metric, summary_stat, sep = '_')) %>%
  filter(radius %in% radii, metric_stat %in% c('raw_sd', 'tri_mean', 'roughness_mean')) %>%
  ggplot(aes(x = value, group = radius, fill = factor(radius))) +
  stat_density(aes(y = ..scaled..), alpha = 0.5, geom = 'polygon', position = 'dodge') +
  facet_grid(metric_stat ~ variable, scales = 'free', labeller = labeller(metric_stat = metric_names, variable = variable_names)) +
  scale_fill_brewer(name = 'Radius (km)', palette = 'Set1') +
  scale_x_continuous(name = 'Metric value') +
  theme_bw() + theme(strip.background = element_rect(fill = NA))

ggsave(file.path(fpfig, 'density_all.png'), p_dens_all, height = 9, width = 10, dpi = 400)



# Do each column separately
x_limits <- list(c(0, 110), c(0, 1400), c(0, 75))

# Maxima for normalizing
fia_melt2 <- melt(fia_vars, id.vars = 1:5, measure.vars = 6:9, variable.name = 'summary_stat') %>% 
  mutate(metric_stat = paste(metric, summary_stat, sep = '_')) %>%
  filter(radius %in% radii, metric_stat %in% c('raw_sd', 'tri_mean', 'roughness_mean'))

fia_maxes <- fia_melt2 %>%
  group_by(variable, metric_stat) %>%
  summarize(maxval = max(value, na.rm = TRUE))

fia_melt2 <- fia_melt2 %>%
  left_join(fia_maxes) %>%
  mutate(scale_value = value/maxval)

p_columns <- list()
for (i in 1:length(geo_vars)) {
  p_columns[[i]] <-
    fia_melt2 %>%
    filter(variable == geo_vars[i])  %>%
    ggplot(aes(x = scale_value, group = radius, fill = factor(radius))) +
    geom_density(aes(y = ..scaled..), alpha = 0.5) +
    facet_grid(metric_stat ~ variable, scales = 'free', labeller = labeller(metric_stat = metric_names, variable = variable_names)) +
    scale_fill_brewer(name = 'Radius (km)', palette = 'Set1') +
    scale_x_continuous(name = 'scaled metric value') +
    scale_y_continuous(name = 'scaled density', expand = c(0,0), limits = c(0,1.05)) +
    theme_bw() + theme(strip.background = element_rect(fill = NA),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
  
  # If i is 2 or 3, remove the axis ticks and labels
  if (i > 1) {
    p_columns[[i]] <- p_columns[[i]] + 
      theme(axis.title.y = element_blank())
  }
  # If i is 1 or 2, remove the y strip and the legend
  if (i < 3) {
    p_columns[[i]] <- p_columns[[i]] +
      theme(strip.text.y = element_blank(),
            legend.position = 'none')
  }
}

# Arrange the three column plots in one figure
p_3col <- cowplot::plot_grid(plotlist = p_columns, rel_widths = c(1.2, 1, 1.85), ncol = 3)
ggsave(file.path(fpfig, 'newsize_density_all_scaled.png'), p_3col, height = 6, width = 6.7, dpi = 400)

# Just elevation in a column
p_elevcol <- fia_melt2 %>%
  filter(variable == 'elevation')  %>%
  ggplot(aes(x = value, group = radius, fill = factor(radius))) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +
  facet_grid(metric_stat ~ variable, scales = 'free', labeller = labeller(metric_stat = metric_names, variable = variable_names)) +
  scale_fill_brewer(name = 'Radius (km)', palette = 'Set1') +
  scale_x_continuous(name = 'metric value', limits = c(0, 1400)) +
  scale_y_continuous(name = 'scaled density', expand = c(0,0), limits = c(0,1.05)) +
  theme_bw() + theme(strip.background = element_rect(fill = NA),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     strip.text.x = element_blank())
ggsave(file.path(fpfig, 'newsize_density_elev.png'), p_elevcol, height = 6, width = 5.5, dpi = 400)

# Demo of circle with plots and geodiversity ------------------------------

# Load FIA data and 5km elevation layer

library(rgdal)
library(raster)

fpras <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/geodata'
ras_elev <- raster(file.path(fpras, 'conus_5k_dem.tif')) # In Albers
ras_tri <- raster(file.path(fpras, 'conus_5k_dem_TRI.tif'))
ras_rough <- raster(file.path(fpras, 'conus_5k_dem_roughness.tif'))

fiafuz <- read.csv('~/FIA/fiafuzzedbyQ.csv')

set.seed(1)
n <- sample(nrow(fiafuz),1)

lonlim <- fiafuz$lonfuzz_aea[n] + c(-30e3, 30e3)
latlim <- fiafuz$latfuzz_aea[n] + c(-30e3, 30e3)

coord_sub <- subset(fiafuz,
                    lonfuzz_aea >= lonlim[1] & lonfuzz_aea <= lonlim[2] & 
                      latfuzz_aea >= latlim[1] & latfuzz_aea <= latlim[2])

anno_circle <- function(xc, yc, r, ...) annotate("path",
                                            x=xc+r*cos(seq(0,2*pi,length.out=1001)),
                                            y=yc+r*sin(seq(0,2*pi,length.out=1001)), ...)

anno_star <- function(xc, yc, ...) annotate("point",
                                            x = xc, y = yc, ...)
 
pts_elev <- rasterToPoints(ras_elev) %>%
  as.data.frame %>%
  filter(between(x, lonlim[1], lonlim[2]), between(y, latlim[1], latlim[2]))
pts_tri <- rasterToPoints(ras_tri) %>%
  as.data.frame %>%
  filter(between(x, lonlim[1], lonlim[2]), between(y, latlim[1], latlim[2]))
pts_rough <- rasterToPoints(ras_rough) %>%
  as.data.frame %>%
  filter(between(x, lonlim[1], lonlim[2]), between(y, latlim[1], latlim[2]))

theme_schem <- theme(legend.position = 'bottom', 
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

ggplot(as.data.frame(pts_elev), aes(x=x, y=y)) +
  geom_tile(aes(fill = conus_5k_dem)) +
  geom_point(data = coord_sub, aes(x=lonfuzz_aea, y=latfuzz_aea), color = 'black') +
  coord_equal(xlim = lonlim + c(5e3,-5e3), ylim = latlim + c(5e3,-5e3)) +
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 20e3, size = 2) +
  anno_star(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], size = 5, pch = 18) +
  scale_fill_gradientn(name = 'Elevation (m)', colours = RColorBrewer::brewer.pal(9,'Oranges')[1:7]) +
  theme_schem

ggsave(file.path(fpfig, 'elevation_demo.png'), height = 3.3, width = 3, dpi = 400) 


p1 <- ggplot(as.data.frame(pts_elev), aes(x=x, y=y)) +
  geom_tile(aes(fill = conus_5k_dem)) +
  coord_equal(xlim = lonlim + c(5e3,-5e3), ylim = latlim + c(5e3,-5e3)) +
  scale_fill_gradientn(name = 'Elevation (m)', colours = RColorBrewer::brewer.pal(9,'Oranges')[1:7]) +
  theme_schem
p2 <- ggplot(as.data.frame(pts_tri), aes(x=x, y=y)) +
  geom_tile(aes(fill = conus_5k_dem_TRI)) +
  coord_equal(xlim = lonlim + c(5e3,-5e3), ylim = latlim + c(5e3,-5e3)) +
  scale_fill_gradientn(name = 'Terrain\nruggedness index', colours = RColorBrewer::brewer.pal(9,'Oranges')[1:7]) +
  theme_schem
p3 <- ggplot(as.data.frame(pts_rough), aes(x=x, y=y)) +
  geom_tile(aes(fill = conus_5k_dem_roughness)) +
  coord_equal(xlim = lonlim + c(5e3,-5e3), ylim = latlim + c(5e3,-5e3)) +
  scale_fill_gradientn(name = 'Roughness', colours = RColorBrewer::brewer.pal(9,'Oranges')[1:7]) +
  theme_schem
ggsave(file.path(fpfig, 'schem_elev.png'), p1, height = 3.3, width = 3, dpi = 400)
ggsave(file.path(fpfig, 'schem_tri.png'), p2, height = 3.3, width = 3, dpi = 400)
ggsave(file.path(fpfig, 'schem_rough.png'), p3, height = 3.3, width = 3, dpi = 400)

p1circ <- p1 + 
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 5e3, size = 1) +
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 10e3, size = 1) +
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 20e3, size = 1) +
  anno_star(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], size = 5, pch = 18)

p2circ <- p2 + 
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 5e3, size = 1) +
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 10e3, size = 1) +
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 20e3, size = 1) +
  anno_star(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], size = 5, pch = 18)

p3circ <- p3 + 
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 5e3, size = 1) +
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 10e3, size = 1) +
  anno_circle(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], 20e3, size = 1) +
  anno_star(fiafuz$lonfuzz_aea[n], fiafuz$latfuzz_aea[n], size = 5, pch = 18)

ggsave(file.path(fpfig, 'schem_elev_circ.png'), p1circ, height = 3.3, width = 3, dpi = 400)
ggsave(file.path(fpfig, 'schem_tri_circ.png'), p2circ, height = 3.3, width = 3, dpi = 400)
ggsave(file.path(fpfig, 'schem_rough_circ.png'), p3circ, height = 3.3, width = 3, dpi = 400)

