# Calculate alpha diversity at all radii for both fia and bbs. 
# I've already calculated the alpha diversity at each plot. 
# Load alpha diversity, for each plot calculate the distance to all the neighbors (split by year if bbs)
# Subset by radius, then calculate the median.
# Edited 19 May: use the new alpha diversity from the HPCC.

fia_alphadiv <- read.csv('/mnt/research/nasabio/data/fia/fia_alphadiv.csv')

radii <- c(5, 10, 20, 50, 100, 150, 200, 300, 400, 500)

library(sp)
library(dplyr)

neighbordiv <- function(x, plot_diversity=fia_alphadiv) {
  neighbordists <- spDistsN1(pts = cbind(plot_diversity$lon, plot_diversity$lat), pt = c(x$lon, x$lat), longlat = TRUE)
  commdat <- list()
  for (i in 1:length(radii)) {
    neighbors_incircle <- plot_diversity[neighbordists <= radii[i], ]
    commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], richness = median(richness,na.rm=T), shannon = median(shannon, na.rm=T), MPD = median(MPD_z, na.rm=T), MNTD = median(MNTD_z, na.rm=T), MPDfunc = median(MPDfunc_z, na.rm=T), MNTDfunc = median(MNTDfunc_z, na.rm=T)))
  }
  as.data.frame(do.call('rbind', commdat))
}

fia_alpha <- fia_alphadiv %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>% 
  do(neighbordiv(.)) # Takes 4 minutes on my machine.

write.csv(fia_alpha, '/mnt/research/nasabio/data/fia/fia_alpha.csv', row.names = FALSE)
