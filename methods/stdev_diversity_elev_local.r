# Run the proof of concept script with lower-resolution elevation data locally.
# Use a single SRTM tile and get some points from it. (90 m resolution)
# Extract different buffer sizes (radii) up to ~100km
# Run different rounding truncations and calculate (1) SD (2) Richness (3) Shannon entropy for each one
# QDR 19 Sep 2017 nasabioxgeo

library(raster)
library(dplyr)

fiapnw <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1]) %>%
  ungroup

# Get a subset of fia coordinates that fall in one srtm tile
fiacoords <- filter(fiacoords, lat > 44.2, lat < 44.8, lon > -122.8, lon < -122.2)
getData(name = 'SRTM', download = TRUE, path = '~/R/srtm/', lon = -122.5, lat = 44.5)

r <- raster('~/R/srtm/srtm_12_04.tif')

summ_stats <- function(dat, truncations) {
  require(vegan)
  dat <- na.omit(dat)
  richnesses <- sapply(truncations, function(x) length(unique(plyr::round_any(dat, x))))
  diversities <- sapply(truncations, function(x) diversity(table(plyr::round_any(dat, x)), index = 'shannon'))
  sds <- sapply(truncations, function(x) sd(plyr::round_any(dat, x)))
  data.frame(truncation = truncations, richness = richnesses, diversity = diversities, sd = sds)
}


radii <- c(5, 10, 20, 50, 100)

all_stats <- list()

pb <- txtProgressBar(0, length(radii), style = 3)

for (i in 1:length(radii)) {
setTxtProgressBar(pb, i)
er <- extract(r, fiacoords[,c('lon','lat')], buffer = radii[i] * 1000)
stat_r <- lapply(er, summ_stats, truncations = c(0.01, 0.1, 1, 10, 100))
all_stats[[length(all_stats) + 1]] <- data.frame(radius = radii[i], lat = rep(fiacoords$lat, each=5), lon=rep(fiacoords$lon, each=5), do.call('rbind', stat_r))
}

close(pb)

# Visualize
library(cowplot)
all_stats <- do.call('rbind', all_stats)

ggplot(all_stats, aes(x = richness, y = sd)) +
  geom_hex() +
  facet_grid(truncation ~ radius) +
  stat_smooth(method = lm, color = 'red', se = FALSE)