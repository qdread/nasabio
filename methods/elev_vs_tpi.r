# Proof of concept: TPI versus sd of elevation
# range of TPI shows the variability in pairwise differences between points and their surrounding elevations. TPI is not absolute value so the mean TPI for a large enough area should be around zero, as peaks should cancel valleys.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/'
#fp <- '/mnt/research/nasabio/data/bbs'

ed <- read.csv(file.path(fp, 'bbs_geodiversity_stats.csv'))

library(dplyr)
library(reshape2)

# Massage so that all summary stats for elevation and TPI each have their own column.
elev <- ed %>%
  filter(variable %in% c('elevation','TPI')) %>%
  select(-richness_geodiv, -diversity_geodiv) %>%
  mutate(range = max - min, cv = mean/sd) %>%
  melt(id.vars = c('rteNo','lon','lat','lon_aea','lat_aea','radius','variable'), variable.name='summary_stat') %>%
  dcast(rteNo + lon + lat + lon_aea + lat_aea + radius ~ variable + summary_stat)

library(ggplot2)

p <- ggplot(elev) +
  facet_wrap(~ radius) +
  theme_bw()

p + geom_point(aes(x = elevation_sd, y = TPI_range))
p + geom_point(aes(x = elevation_sd, y = TPI_sd))
p + geom_point(aes(x = elevation_sd, y = TPI_max))
p + geom_point(aes(x = elevation_range, y = TPI_range))

# None of these will probably work. TPI_range is too noisy. Some quantile would be better, but I didn't pull quantiles.
# Mean absolute difference is now added as of 20 Nov.
hist(elev$TPI_mean)

p + geom_point(aes(x = elevation_sd, y = TPI_mean))


# Same metrics at different radii -----------------------------------------


# Compare difference in elevation sd across radii, with difference in TPI mean across radii.

# For example, compare elevation sd at 50 km with elevation sd at 200 km.

# Use FIA for now although we do not have the correct TPI mean, because we have the small radii.
# Right now, we just have bioclim and elevation standard deviations.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

ed <- read.csv(file.path(fp, 'fia_elev_stats_unfuzzed.csv'))

library(dplyr)
library(reshape2)
library(ggplot2)
library(GGally)

# Put in wide format with different columns for different radii

elev_by_rad <- ed %>%
  select(-variable) %>%
  mutate(range = max - min, cv = mean/sd, radius = paste('r', radius, sep = '_')) %>%
  filter(!is.na(mean)) %>%
  melt(id.vars = c('PLT_CN', 'STATECD', 'COUNTYCD', 'PLOT', 'radius'), variable.name='summary_stat') %>%
  dcast(PLT_CN + STATECD + COUNTYCD + PLOT ~ radius + summary_stat)

sds_by_rad <- elev_by_rad %>% 
  select(contains('sd'))

# Plot paired plots
my_hex <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_hex(...) + 
    scale_fill_continuous(low = 'gray90', high = 'black')
  p
}


pdf('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots/elevsd_radius_pairplot.pdf', height = 10, width = 10)
ggpairs(sds_by_rad[, c('r_5_sd', 'r_10_sd', 'r_20_sd', 'r_50_sd', 'r_100_sd')],
        diag = list(continuous = wrap('barDiag', bins=20)),
        lower = list(continuous = my_hex)) + 
  theme_bw()
dev.off()