# Split up bird data by longevity
# QDR / NASAbioXgeo / 28 Jun 2018

google_drive <- 'C:/Users/Q/google_drive/NASABiodiversityWG' # change as needed

library(dplyr)
library(ggplot2)

# Load full imputed dataset
traits <- read.csv(file.path(google_drive, 'birdXtree_ScalingPaper/bird_trait_data/bbs_traits_imputed.csv'), stringsAsFactors = FALSE)

# Look at longevity numbers with histograms/density plots

ggplot(traits, aes(x = longevity_y, y = maximum_longevity_y)) + geom_point() + theme_minimal() # Quite a bit of variation

# Correct the maximum longevities that are less than the typicals
traits <- traits %>%
  mutate(maximum_longevity_y = pmax(maximum_longevity_y, longevity_y))

# Typical
ggplot(traits, aes(x = longevity_y)) + geom_histogram() + theme_minimal() # Regular scale
ggplot(traits, aes(x = longevity_y)) + geom_histogram() + theme_minimal() + scale_x_log10(breaks = c(1,3,10,30)) # Log scale

# Maximum
ggplot(traits, aes(x = maximum_longevity_y)) + geom_histogram() + theme_minimal() # Regular scale
ggplot(traits, aes(x = maximum_longevity_y)) + geom_histogram() + theme_minimal() + scale_x_log10(breaks = c(1,3,10,30)) # Log scale

# There might be a slight bimodality of maximum longevity, when log transformed.

longevity_cluster <- kmeans(na.omit(log10(traits$maximum_longevity_y)), centers = 2)

traits$longevity_group <- NA
traits$longevity_group[!is.na(traits$maximum_longevity_y)] <- longevity_cluster$cluster

min(traits$maximum_longevity_y[traits$longevity_group ==2], na.rm=TRUE) # 14.5 years is split

plot_longevity <- traits %>%
  mutate(longevity_group = factor(longevity_group)) %>%
  ggplot(aes(x = maximum_longevity_y, group = longevity_group, fill = longevity_group)) + 
  geom_histogram() + theme_minimal() + scale_x_log10(breaks = c(1,3,10,30)) +
  scale_fill_discrete(name = 'Maximum longevity', labels = c('<= 14.5 years', '> 14.5 years')) +
  theme(legend.position = 'bottom')

# Try time to maturity
# Plot female vs male time to maturity to see if it matters which is used
ggplot(traits, aes(x = female_maturity_d, y = male_maturity_d)) + geom_point() + theme_minimal() # Very similar, should not matter which is used

ggplot(traits, aes(x = female_maturity_d)) + geom_histogram() + theme_minimal() # Regular scale
ggplot(traits, aes(x = female_maturity_d)) + geom_histogram() + theme_minimal() + scale_x_log10(breaks = c(100,365,1000)) # Log scale

table(traits$female_maturity_d <= 365) # Actually not a bad split.

plot_maturity <- ggplot(traits, aes(x = female_maturity_d, group = female_maturity_d <= 365, fill = female_maturity_d <= 365)) + 
  geom_histogram() + 
  theme_minimal() + 
  scale_x_log10(breaks = c(100,365,1000)) +
  scale_fill_discrete(name = 'Female maturity', labels = c('<= 1 year', '> 1 year')) +
  theme(legend.position = 'bottom')

pdf(file.path(google_drive, 'birdXtree_ScalingPaper/bbs_groups.pdf'), height=6, width=6)
plot_longevity
plot_maturity
dev.off()

with(traits, table(longevity_group, female_maturity_d <= 365))
