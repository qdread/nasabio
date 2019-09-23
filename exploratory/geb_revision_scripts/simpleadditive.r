# Simple script to use the very basic definition of additive beta diversity as gamma - mean(alpha)

fp <- '~/Dropbox/projects/nasabiodiv'

bbsa <- read.csv(file.path(fp, 'bbs_alpha_1year.csv'))
bbsb <- read.csv(file.path(fp, 'bbs_betatdpdfd_1year.csv'))
bbsg <- read.csv(file.path(fp, 'bbs_gamma_1year.csv'))

library(tidyverse)

# Get the very basic additive beta diversity.
bbsadditive <- bbsa %>%
  select(rteNo, lon, lat, radius, richness) %>%
  rename(richness_alpha = richness) %>%
  left_join(bbsg %>%
              select(rteNo, lon, lat, radius, richness) %>%
              rename(richness_gamma = richness)) %>%
  mutate(richness_beta = richness_gamma - richness_alpha) %>%
  filter(radius == 50)

# Join the sorensen beta to this.
bbsadditive <- bbsadditive %>% left_join(bbsb %>% filter(radius == 50) %>% select(rteNo, lon, lat, radius, beta_td_sorensen_pa))

library(GGally)
ggpairs(bbsadditive[, 5:8])

# Do the same with FIA.
load(file.path(fp, 'modelfits/fia_spatial_mm_dat_50k.RData'))

fiaadditive <- fiabio %>% 
  select(PLT_CN, alpha_richness, beta_td_sorensen_pa, gamma_richness) %>%
  mutate(beta_richness = gamma_richness - alpha_richness)

ggpairs(fiaadditive[,-1])
