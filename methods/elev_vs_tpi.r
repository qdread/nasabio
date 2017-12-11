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

# Edited 11 Dec to compare elevation with the absolute topographic position mean.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

ed <- read.csv(file.path(fp, 'fia_pnw_elev_stats.csv'))

library(dplyr)
library(reshape2)
library(ggplot2)
library(GGally)

# Put in wide format with different columns for different radii

elev_by_rad <- ed %>%
  mutate(range = max - min, cv = mean/sd, radius = paste('r', radius, sep = '_')) %>%
  filter(!is.na(mean), variable == 'elevation') %>%
  select(-variable) %>%
  melt(id.vars = c('PLT_CN', 'STATECD', 'COUNTYCD', 'PLOT', 'radius'), variable.name='summary_stat') %>%
  dcast(PLT_CN + STATECD + COUNTYCD + PLOT ~ radius + summary_stat)

sds_by_rad <- elev_by_rad %>% 
  select(contains('sd'))

tpi_by_rad <- ed %>%
  mutate(range = max - min, cv = mean/sd, radius = paste('r', radius, sep = '_')) %>%
  filter(!is.na(mean), variable == 'TPI') %>%
  select(-variable) %>%
  melt(id.vars = c('PLT_CN', 'STATECD', 'COUNTYCD', 'PLOT', 'radius'), variable.name='summary_stat') %>%
  dcast(PLT_CN + STATECD + COUNTYCD + PLOT ~ radius + summary_stat)

tpimean_by_rad <- tpi_by_rad %>%
  select(contains('mean'))

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

# Comparison of SD elevation and mean TPI at different radii
elev_tpi_plots <- list()
for (i in c(5,10,20,50,100)) {
  dat <- data.frame(elev_sd = sds_by_rad[,paste('r',i,'sd',sep='_')],
                    tpi_mean = tpimean_by_rad[,paste('r',i,'mean',sep='_')])
  r2 <- summary(lm(tpi_mean~elev_sd, data=dat))$r.sq
  plot_i <- ggplot(dat, aes(x = elev_sd, y = tpi_mean)) +
    geom_hex() +
    geom_text(data = data.frame(elev_sd=-Inf, tpi_mean=Inf, lab=paste('R^2 ==', round(r2,2))), aes(label=lab), parse = TRUE, hjust = -1, vjust = 1) +
    scale_fill_continuous(low = 'gray90', high = 'black') +
    theme_bw() +
    ggtitle(paste(i, 'km')) +
    labs(x = 'Elevation SD', y = 'Absolute topographic difference mean') +
    theme(legend.position = 'none')
  elev_tpi_plots[[length(elev_tpi_plots) + 1]] <- plot_i
}

library(gridExtra)
pdf('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots/tpi_mean_vs_elev_sd.pdf', height = 4, width = 12)
grid.arrange(grobs = elev_tpi_plots, nrow = 1)
dev.off()