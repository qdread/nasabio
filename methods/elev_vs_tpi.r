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
# What we need is the mean ABSOLUTE difference, because mean TPI is ~ 0.
hist(elev$TPI_mean)