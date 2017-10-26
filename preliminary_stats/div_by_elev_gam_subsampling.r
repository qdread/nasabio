# Conceptual paper model fitting: FIA a,b,g diversity ~ elevation SD

# Procedure in pseudocode:
# Load a,b,g diversity data and elevation diversity data
# Throw out all points within 100 km of Mexico and Canada.
# for radius in 5, 10, 20, 50, 100 {
#   for iteration in 1:999 {
#     for y_variable in alpha, beta, gamma {
#       Do iterative search to get ~20 points from the dataset (since that is as many 100 km radius ones as we can get)
#       Fit GAM to those points
#       Save result to some sort of array
# }}}


# NOTE: All needed CSVs are on the hpcc but I have downloaded them locally because it is faster. 
# Change file path to the second file path to get files from hpcc.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_elev_stats_unfuzzed.csv'))
ad <- read.csv(file.path(fp, 'fia_alpha.csv'))
gd <- read.csv(file.path(fp, 'fia_gammadiv.csv'))

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_richness','beta_richness','gamma_richness')

# Combine into a single data frame.
library(dplyr)
biogeo <- ed %>%
  select(PLT_CN, radius, sd) %>%
  filter(radius %in% radii) %>%
  rename(elevation_sd = sd) %>%
  left_join(ad %>% select(PLT_CN, radius, richness) %>% rename(alpha_richness = richness)) %>%
  left_join(gd %>% select(PLT_CN, radius, richness) %>% rename(gamma_richness = richness))

# Add latitude and longitudes from unfuzzed (on local drive only)
fiacoords <- read.csv('~/FIA/pnw.csv')

biogeo <- biogeo %>%
  left_join(fiacoords %>% select(CN, ACTUAL_LAT, ACTUAL_LON) %>% rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON))

# Function for iterative search.
source('~/GitHub/nasabio/methods/SRS_iterative.r')

n_iter <- 999

# R2 of gam. array with number of diversity types (a,b,g) by number of radii by number of iterations
r2_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))

for (i in 1:length(div_names)) {
  for (j in 1:length(radii)) {
    for (k in 1:n_iter) {
      # Subsample
      
      # Fit GAM
      
      # Save r-squared
    }
  }
}