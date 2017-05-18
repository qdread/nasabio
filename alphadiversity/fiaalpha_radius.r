# Calculate alpha diversity at all radii for both fia and bbs. 
# I've already calculated the alpha diversity at each plot. 
# Load alpha diversity, for each plot calculate the distance to all the neighbors (split by year if bbs)
# Subset by radius, then calculate the median.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

load(file.path(fp, 'fia_diversitymetrics.RData'))

# Also add functional diversity to this.
fia_fd <- read.csv(file.path(fp, 'fia_fd.csv'))
plot_diversity <- cbind(plot_diversity, fia_fd)

# Calculate plot diversity within radii

library(sp)
library(dplyr)

fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = F)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

plot_diversity <- left_join(plot_diversity, fiacoords)
radii <- c(5, 10, 20, 50, 100, 150, 200, 300, 400, 500)


neighbordiv <- function(x) {
  neighbordists <- spDistsN1(pts = cbind(plot_diversity$lon, plot_diversity$lat), pt = c(x$lon, x$lat), longlat = TRUE)
  commdat <- list()
  for (i in 1:length(radii)) {
    neighbors_incircle <- plot_diversity[neighbordists <= radii[i], ]
    commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], richness = median(richness,na.rm=T), shannon_basalarea = median(shannon_basalarea, na.rm=T), mpd = median(mpd_z, na.rm=T), mntd = median(mntd_z, na.rm=T)))
  }
  as.data.frame(do.call('rbind', commdat))
}

fia_alpha <- plot_diversity %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>% 
  do(neighbordiv(.)) # Takes 4 minutes on my machine.

write.csv(fia_alpha, file.path(fp, 'fia_alpha.csv'), row.names = FALSE)
