# Calculate median alpha-diversity of bbs by radius.

# Load bbs alpha diversity and route coordinates
load('/mnt/research/aquaxterra/DATA/raw_data/BBS/bbs_div_byroute_presence.r')

# Edited 05 May: load the correct coordinates and replace the wrong ones.
bbsll <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv')

library(dplyr)

bbs_div_byroute <- bbs_div_byroute %>% 
	select(-latitude, -longitude) %>%
	left_join(bbsll)

# For each year and route number, get the median alpha diversity within each radius.
radii <- c(50, 75, 100) # in km

library(sp)

bbs_div_byroute_filtered <- bbs_div_byroute %>% filter(!is.na(lat), !is.na(lon), year >= 1997)

neighbordiv <- function(x) {
	neighborrows <- bbs_div_byroute_filtered[bbs_div_byroute_filtered$year == x$year, ]
	neighbordists <- spDistsN1(pts = cbind(neighborrows$lon, neighborrows$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- neighborrows[neighbordists <= radii[i], ]
		commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], richness = median(richness), FRic = median(FRic), FEve = median(FEve), FDiv = median(FDiv), FDis = median(FDis), PD = median(PD), mpd.obs.z = median(mpd.obs.z), mntd.obs.z = median(mpd.obs.z)))
	}
	as.data.frame(do.call('rbind', commdat))
}

bbs_alpha <- bbs_div_byroute_filtered %>% 
	group_by(year, rteNo) %>%
	do(neighbordiv(.))
	
bbs_alpha <- bbs_alpha %>% ungroup %>% left_join(bbs_div_byroute_filtered %>% select(year,rteNo,lat,lon))
	
write.csv(bbs_alpha, file = '/mnt/research/nasabio/data/bbs/bbs_alpha.csv', row.names = FALSE)	
	