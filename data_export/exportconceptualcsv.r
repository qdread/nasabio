# Put together CSV of all data used in conceptual paper, with fuzzed locations
# Includes only Pacific Northwest

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_elev_stats_unfuzzed.csv'))
ad <- read.csv(file.path(fp, 'fia_alpha.csv'))
bd <- read.csv(file.path(fp, 'fia_beta.csv'))
gd <- read.csv(file.path(fp, 'fia_gammadiv.csv'))

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_diversity','beta_diversity','gamma_diversity')

library(dplyr)

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

# Get rid of Alaska
biogeo <- biogeo %>%
  filter(STATECD != 2)

# Join with fuzzed coordinates
fiaraw <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = F) 
fuzzedcoords <- fiaraw %>%
  group_by(PLT_CN, STATECD, COUNTYCD, PLOT) %>%
  summarize(lat_fuzz = LAT_FUZZSWAP[1], lon_fuzz = LON_FUZZSWAP[1])
biogeo <- biogeo %>%
  left_join(fuzzedcoords) %>%
  dplyr::select(PLT_CN, STATECD, COUNTYCD, PLOT, lat_fuzz, lon_fuzz, everything())

# Retain only diversity columns and exponentiate the alpha and gamma ones to get effective species number
biogeo <- biogeo %>%
  dplyr::select(-contains('richness')) %>%
  mutate(alpha_diversity = exp(alpha_diversity),
         gamma_diversity = exp(gamma_diversity))

biogeo <- biogeo %>%
  arrange(PLT_CN, radius)

write.csv(biogeo, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/conceptual_paper_alldata.csv', row.names = FALSE)
