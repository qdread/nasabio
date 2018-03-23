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
  filter(grepl('bio1_5k', variable))

fia_vars <- rbind(fia_elev_long, fia_other, fia_bio5k) 
write.csv(fia_vars, file = '/mnt/research/nasabio/temp/fia_vars_iale.csv', row.names = FALSE)

fpbbs <- '/mnt/research/nasabio/data/bbs'

bbs_vars <- read.csv(file.path(fpbbs, 'bbs_geodiversity.csv'), stringsAsFactors = FALSE)

# Keep only elevation, human footprint, and mean annual precipitation.

bbs_vars <- bbs_vars %>%
  filter(grepl('elevation_5k|human_footprint_5k|bio1_5k', variable))
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
  select(-richness_geodiv, -diversity_geodiv, -mode, -lon, -lat, -lon_aea, -lat_aea)

bbs_vars_melt <- melt(bbs_vars, id.vars = 1:5, variable.name = 'summary_stat', value.name = 'value')
bbs_vars_cast <- dcast(bbs_vars_melt, rteNo + variable + resolution + radius ~ metric + summary_stat)
bbs_vars_cast_wide <- dcast(bbs_vars_melt, rteNo + variable + resolution ~ metric + summary_stat + radius)

# FIA (a lot bigger, csv is 1gb)
fia_vars <- read.csv(file.path(fp, 'fia_vars_iale.csv'), stringsAsFactors = FALSE)

fia_vars <- fia_vars %>%
  mutate(variable = gsub('human_','human',variable)) %>%
  separate(variable, into = c('variable','resolution','metric')) %>%
  mutate(metric = if_else(is.na(metric), 'raw', metric)) 

fia_vars_melt <- melt(fia_vars, id.vars = 1:5, variable.name = 'summary_stat', value.name = 'value')
fia_vars_cast <- dcast(fia_vars_melt, PLT_CN + variable + resolution + radius ~ metric + summary_stat)
fia_vars_cast_wide <- dcast(fia_vars_melt, PLT_CN + variable + resolution ~ metric + summary_stat + radius)

rm(fia_vars, fia_vars_melt)

# Make plots comparing the metrics at different radii ---------------------

library(ggplot2)

radii <- c(5, 10, 20, 50, 100, 200)

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  ggplot(aes(x = raw_sd, y = roughness_mean)) +
    facet_wrap(~ radius) +
    geom_hex() +
    scale_fill_gradient(low = 'gray80', high = 'black') +
    theme_bw() +
    theme(strip.background = element_blank()) +
    ggtitle('Roughness vs. standard deviation', 'Elevation') +
    labs(x = 'Standard deviation', y = 'Roughness')
 
fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  ggplot(aes(x = raw_sd, y = tri_mean)) +
    facet_wrap(~ radius) +
    geom_hex() +
    scale_fill_gradient(low = 'gray80', high = 'black') +
    theme_bw() +
    theme(strip.background = element_blank()) +
    ggtitle('Terrain ruggedness index vs. standard deviation', 'Elevation') +
    labs(x = 'Standard deviation', y = 'TRI')   

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  ggplot(aes(x = roughness_mean, y = tri_mean)) +
    facet_wrap(~ radius) +
    geom_hex() +
    scale_fill_gradient(low = 'gray80', high = 'black') +
    theme_bw() +
    theme(strip.background = element_blank()) +
    ggtitle('Terrain ruggedness index vs. roughness', 'Elevation') +
    labs(x = 'Roughness', y = 'TRI')   


# Make plot showing how a single metric changes with radius ---------------

set.seed(555)
plot_sample <- sample(unique(fia_vars_cast_wide$PLT_CN), 1000)

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

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  ggplot(aes(x = raw_sd, group = radius, fill = factor(radius))) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(name = 'Standard deviation') +
    scale_fill_brewer(palette = 'Set1') +
    theme_bw() +
    ggtitle('Standard deviation change with radius', 'Elevation')

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  ggplot(aes(x = roughness_mean, group = radius, fill = factor(radius))) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(name = 'Roughness') +
    scale_fill_brewer(palette = 'Set1') +
    theme_bw() +
    ggtitle('Roughness change with radius', 'Elevation')

fia_vars_cast %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  ggplot(aes(x = tri_mean, group = radius, fill = factor(radius))) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(name = 'TRI') +
    scale_fill_brewer(palette = 'Set1') +
    theme_bw() +
    ggtitle('Terrain ruggedness index change with radius', 'Elevation')
