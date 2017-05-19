# Calculate median alpha-diversity of bbs by radius.
# edit 19 may: new correct diversity.

# Load bbs alpha diversity and route coordinates
bbs_alphadiv <- read.csv('/mnt/research/nasabio/data/bbs/bbs_alphadiv.csv')

library(dplyr)

# For each year and route number, get the median alpha diversity within each radius.
radii <- c(50, 75, 100, 150, 200, 300, 400, 500) # in km

library(sp)

bbs_alphadiv_filtered <- bbs_alphadiv %>% filter(!is.na(lat), !is.na(lon), year >= 1997)

neighbordiv <- function(x) {
	neighborrows <- bbs_alphadiv_filtered[bbs_alphadiv_filtered$year == x$year, ]
	neighbordists <- spDistsN1(pts = cbind(neighborrows$lon, neighborrows$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- neighborrows[neighbordists <= radii[i], ]
		commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], richness = median(richness), MPD_pa_z = median(MPD_pa_z, na.rm=TRUE), MNTD_pa_z = median(MNTD_pa_z, na.rm=TRUE), MPDfunc_pa_z = median(MPDfunc_pa_z, na.rm = TRUE), MNTDfunc_pa_z = median(MNTDfunc_pa_z[is.finite(MNTDfunc_pa_z)])))
	}
	as.data.frame(do.call('rbind', commdat))
}

bbs_alpha <- bbs_alphadiv_filtered %>% 
	group_by(year, rteNo) %>%
	do(neighbordiv(.))
	
bbs_alpha <- bbs_alpha %>% ungroup %>% left_join(bbs_alphadiv_filtered %>% select(year,rteNo,lat,lon))
	
write.csv(bbs_alpha, file = '/mnt/research/nasabio/data/bbs/bbs_alpha.csv', row.names = FALSE)	
	