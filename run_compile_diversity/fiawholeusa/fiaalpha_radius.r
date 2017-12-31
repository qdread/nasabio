# Calculate alpha diversity at all radii for both fia and bbs. 
# I've already calculated the alpha diversity at each plot. 
# Load alpha diversity, for each plot calculate the distance to all the neighbors (split by year if bbs)
# Subset by radius, then calculate the median.
# Edited 19 May: use the new alpha diversity from the HPCC.
# Edited 25 October: Use coordinates of new unfuzzed plots.
# Edited 31 December: Whole USA unfuzzed.

alphadiv_list <- list()
for (i in 1:1500) {
	load(paste0('/mnt/research/nasabio/data/fia/diversity/usa/alpha_', i, '.r'))
	alphadiv_list[[i]] <- alpha_div
}

fia_alphadiv <- do.call('rbind', alphadiv_list)
names(fia_alphadiv)[1] <- 'PLT_CN'

write.csv(fia_alphadiv, file = '/mnt/research/nasabio/data/fia/fiausa_alphadiv.csv', row.names = FALSE)

fia_alphadiv <- read.csv('/mnt/research/nasabio/data/fia/fiausa_alphadiv.csv')

load('~/data/fiaworkspace_spatial_wholeusa.r')

radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300)

library(sp)
library(dplyr)

neighbordiv <- function(x, plot_diversity=fia_alphadiv, plot_metadata = fiacoords) {
  neighbordists <- spDistsN1(pts = cbind(plot_metadata$lon, plot_metadata$lat), pt = c(x$lon, x$lat), longlat = TRUE)
  commdat <- list()
  for (i in 1:length(radii)) {
    neighbors_incircle <- plot_diversity[neighbordists <= radii[i], ]
    commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], 
											   richness = median(richness,na.rm=T), 
											   shannon = median(shannon, na.rm=T), 
											   evenness = median(evenness, na.rm=T),
											   MPD = median(MPD_z, na.rm=T), 
											   MNTD = median(MNTD_z, na.rm=T), 
											   MPDfunc = median(MPDfunc_z, na.rm=T), 
											   MNTDfunc = median(MNTDfunc_z, na.rm=T),
											   MPD_pa = median(MPD_pa_z, na.rm=T), 
											   MNTD_pa = median(MNTD_pa_z, na.rm=T), 
											   MPDfunc_pa = median(MPDfunc_pa_z, na.rm=T), 
											   MNTDfunc_pa = median(MNTDfunc_pa_z, na.rm=T)))
  }
  as.data.frame(do.call('rbind', commdat))
}

fia_alpha <- fiacoords %>%
  rowwise %>% 
  do(neighbordiv(.)) # Takes a few hours for the whole USA
  
fia_alpha <- cbind(fiacoords[rep(1:nrow(fiacoords), each = length(radii)), 1:4], fia_alpha)  
  
write.csv(fia_alpha, '/mnt/research/nasabio/data/fia/fiausa_alpha.csv', row.names = FALSE)
