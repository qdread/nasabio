# Calculate alpha diversity at all radii for both fia and bbs. 
# I've already calculated the alpha diversity at each plot. 
# Load alpha diversity, for each plot calculate the distance to all the neighbors (split by year if bbs)
# Subset by radius, then calculate the median.
# Edited 19 May: use the new alpha diversity from the HPCC.
# Edited 25 October: Use coordinates of new unfuzzed plots.
# Edited 31 December: Whole USA unfuzzed.
# Edited 1 March 2018: get the subset of only "natural" plots.
# Edited 28 Nov 2018: update all file paths for updated FIA diversity metrics with macroplots removed.

alphadiv_list <- list()
for (i in 1:1000) {
	load(paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/alpha_', i, '.r'))
	alphadiv_list[[i]] <- alpha_div
}

fia_alphadiv <- do.call('rbind', alphadiv_list)
names(fia_alphadiv)[1] <- 'PLT_CN'

write.csv(fia_alphadiv, file = '/mnt/research/nasabio/data/fia/biodiversity_CSVs/updated_nov2018/fiausa_alphadiv.csv', row.names = FALSE)

fia_alphadiv <- read.csv('/mnt/research/nasabio/data/fia/biodiversity_CSVs/updated_nov2018/fiausa_alphadiv.csv')
plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
fia_alphadiv <- subset(fia_alphadiv, PLT_CN %in% plantation$PLT_CN[!plantation$plantation])

library(dplyr)

load('/mnt/home/qdr/data/fiaworkspace_spatial_wholeusa_2018.r')

fiacoords <- fiacoords %>% left_join(plantation) %>% filter(!plantation)

radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300)
radii <- c(5, 10, 20, 50, 75, 100) # Edited 13 Dec 2018 because new version does not have to go above 100 km

library(sp)

neighbordiv <- function(x, plot_diversity=fia_alphadiv, plot_metadata = fiacoords) {
  neighbordists <- spDistsN1(pts = cbind(plot_metadata$lon, plot_metadata$lat), pt = c(x$lon, x$lat), longlat = TRUE)
  commdat <- list()
  for (i in 1:length(radii)) {
    neighbors_incircle <- plot_diversity[neighbordists <= radii[i], ]
    commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], 
											   richness = median(richness,na.rm=T), 
											   shannon = median(shannon, na.rm=T), 
											   evenness = median(evenness, na.rm=T),
											   MPD_z = median(MPD_z, na.rm=T), 
											   MNTD_z = median(MNTD_z, na.rm=T), 
											   MPDfunc_z = median(MPDfunc_z, na.rm=T), 
											   MNTDfunc_z = median(MNTDfunc_z, na.rm=T),
											   MPD_pa_z = median(MPD_pa_z, na.rm=T), 
											   MNTD_pa_z = median(MNTD_pa_z, na.rm=T), 
											   MPDfunc_pa_z = median(MPDfunc_pa_z, na.rm=T), 
											   MNTDfunc_pa_z = median(MNTDfunc_pa_z, na.rm=T)))
  }
  as.data.frame(do.call('rbind', commdat))
}

fia_alpha <- fiacoords %>%
  rowwise %>% 
  do(neighbordiv(.)) # Takes a few hours for the whole USA
  
fia_alpha <- cbind(PLT_CN = fiacoords[rep(1:nrow(fiacoords), each = length(radii)), 1, drop = FALSE], fia_alpha)  # Don't include the lat and long.
  
write.csv(fia_alpha, '/mnt/research/nasabio/data/fia/biodiversity_CSVs/updated_nov2018/fiausa_natural_alpha.csv', row.names = FALSE)
