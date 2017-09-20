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
  sds <- sd(dat)
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
write.csv(all_stats, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/proofofconceptstats.csv', row.names = FALSE)

ggplot(all_stats, aes(x = sd, y = richness)) +
  geom_point() +
  facet_grid(truncation ~ radius) +
  stat_smooth(method = lm, color = 'red', se = FALSE) +
  panel_border(colour = 'black')

ggplot(all_stats, aes(x = sd, y = diversity)) +
  geom_point() +
  facet_grid(truncation ~ radius) +
  stat_smooth(method = lm, color = 'red', se = FALSE) +
  panel_border(colour = 'black')

# Calculate R^2
all_r2 <- all_stats %>%
  group_by(radius, truncation) %>%
  summarize(r2_richness = summary(lm(richness ~ sd))$r.sq,
            r2_diversity = summary(lm(richness ~ diversity))$r.sq)

# Observations
# 1. The absolute value of richness and diversity are very dependent on the arbitrary threshold of rounding that we choose in order to force a continuous variable to act like a discrete variable
# 2. (Trivially) as the rounding threshold becomes more coarse, you end up with basically no variation in richness or diversity, regardless of the true variation in the elevations.
# 3. The relationship between sd and richness, and between sd and diversity, is different at different radii. I believe this is an artifact of the number of data points in each radius.


# Histogram rules applied to one vector of elevation data
# Return number of bins
sturges_rule <- function(x) {
  k <- ceiling(log(length(x), base = 2)) + 1
  return(k)
}
scott_rule <- function(x) {
  h <- 3.5 * sd(x) * length(x)^(-1/3)
  k <- ceiling(diff(range(x))/h)
  return(k)
}
freedman_rule <- function(x) {
  h <- 2 * IQR(x) * length(x)^(-1/3) # Width of each bin
  k <- ceiling(diff(range(x))/h) # Number of bins of width h that you can split the data into
  return(k)
}

histbinrules <- sapply(1:6, function(x) {
  dat <- sample(er[[1]], 10^x)
  c(sturges = sturges_rule(dat), scott = scott_rule(dat), freedman = freedman_rule(dat))
})

# Redo extractions with Freedman rule

summ_stats_freedman <- function(dat) {
  require(vegan)
  dat <- na.omit(dat)
  n_breaks <- freedman_rule(dat)
  dat_table <- as.numeric(table(cut(dat, breaks = n_breaks)))
  
  data.frame(richness = length(dat_table), diversity = diversity(dat_table, index = 'shannon'), sd = sd(dat))
}


radii <- c(5, 10, 20, 50, 100)

all_stats_freed <- list()

pb <- txtProgressBar(0, length(radii), style = 3)

for (i in 1:length(radii)) {
  setTxtProgressBar(pb, i)
  er <- extract(r, fiacoords[,c('lon','lat')], buffer = radii[i] * 1000)
  stat_r <- lapply(er, summ_stats_freedman)
  all_stats_freed[[length(all_stats_freed) + 1]] <- data.frame(radius = radii[i], lat = fiacoords$lat, lon = fiacoords$lon, do.call('rbind', stat_r))
}

close(pb)


all_stats_freed <- do.call('rbind', all_stats_freed)
write.csv(all_stats_freed, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/proofofconceptstatsfreed.csv', row.names = FALSE)

all_stats_freed <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/code/proofofconceptstatsfreed.csv', stringsAsFactors = FALSE)

library(cowplot)

ggplot(all_stats_freed, aes(x = sd, y = richness)) +
  geom_point() +
  facet_wrap(~ radius) +
  stat_smooth(method = lm, color = 'red', se = FALSE) +
  panel_border(colour = 'black')

ggplot(all_stats_freed, aes(x = sd, y = diversity)) +
  geom_point() +
  facet_wrap(~ radius) +
  stat_smooth(method = lm, color = 'red', se = FALSE) +
  panel_border(colour = 'black')

ggplot(all_stats_freed, aes(x = richness, y = diversity)) +
  geom_point() +
  facet_wrap(~ radius) +
  stat_smooth(method = lm, color = 'red', se = FALSE) +
  panel_border(colour = 'black')
