# Relationships among geodiversity variables at different scales and using different metrics

# Metric: standard deviation, TRI, roughness (mean)
# Scale: small to big

# Use a small subset at first

path_to_file <- 'C:/Users/Q/google_drive/NASABiodiversityWG/SampleData'
##############################################################

fia_dat <- read.csv(file.path(path_to_file, 'fia_subset_wideform.csv'), stringsAsFactors = FALSE)
bbs_dat <- read.csv(file.path(path_to_file, 'bbs_subset_wideform.csv'), stringsAsFactors = FALSE)
metadata <- read.csv(file.path(path_to_file, 'metadata.csv'), stringsAsFactors = FALSE)
source(file.path(path_to_file, 'predictorsubset.r'))

fia_topo_geodiv <- predictor_subset(fia_dat, metadata,
                                          variable_type = 'geo',
                                          radius = 'all',
                                          summary_statistic = c('sd','tri','roughness'),
                                          variable_category = c('topography'),
                                          resolution = 'native')

library(reshape2)
library(ggplot2)
library(GGally)
library(tidyr)
library(dplyr)
fia_topo_long <- melt(fia_topo_geodiv, id.vars = 1:6) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(stat = sapply(strsplit(variable,'_'), function(x) x[length(x)]),
         radius = sapply(strsplit(variable,'_'), function(x) x[length(x)-1]),
         variable = sapply(strsplit(variable,'_'), function(x) paste(x[1:2], collapse = '_')))



fia_elev_sdtrirough <- fia_topo_geodiv %>%
  select(PLT_CN, contains('elevation')) %>%
  select(PLT_CN, contains('sd'), contains('tri'), contains('roughness'), -contains('point')) %>%
  melt(id.vars = 'PLT_CN') %>%
  mutate(variable = as.character(variable)) %>%
  mutate(radius = sapply(strsplit(variable,'_'), function(x) x[length(x)-1]),
         variable = ifelse(grepl('tri', variable), 'tri',ifelse(grepl('rough',variable), 'roughness', 'sd'))) %>%
  dcast(PLT_CN + radius ~ variable)
         
# TRI vs Std Dev  
ggplot(fia_elev_sdtrirough, aes(x = sd, y = tri)) + 
  facet_wrap(~ as.numeric(radius)) +
  geom_point() +
  theme_bw()

# TRI vs Roughness
ggplot(fia_elev_sdtrirough, aes(x = roughness, y = tri)) + 
  facet_wrap(~ as.numeric(radius)) +
  geom_point() +
  theme_bw()

# TRI at 5 km vs TRI at 100 km for all the points
fia_elev_sdtrirough %>%
  filter(radius %in% c(5,100)) %>%
  dcast(PLT_CN ~ radius, value.var = 'tri') %>%
  setNames(c('PLT_CN','tri_100','tri_5')) %>%
  ggplot(aes(x = tri_5, y = tri_100)) + geom_point() + theme_bw()

fia_elev_sdtrirough %>%
  filter(radius %in% c(5,100)) %>%
  dcast(PLT_CN ~ radius, value.var = 'sd') %>%
  setNames(c('PLT_CN','sd_100','sd_5')) %>%
  ggplot(aes(x = sd_5, y = sd_100)) + geom_point() + theme_bw()