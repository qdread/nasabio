# Compile alpha diversity by radius: BBS within route
# ---------------------------------------------------

# Modified 28 June: include 1 kilometer (middle two stops)

# Get alpha averages for the three stop subsets within each route, for BBS

bbs_alphadiv <- read.csv('/mnt/research/nasabio/data/bbs/biodiversity_CSVs/withinroute/bbs_withinroute_alphabystop.csv', stringsAsFactors = FALSE)

library(dplyr)

# For each year and route number, get the median alpha diversity within each radius.
radii <- c(1, 5, 10, 20) # in km
stop_bounds <- rbind(c(25,26), c(20,31), c(13,37), c(1,50)) # First and last stops to be used for each radius (approx. 5,10,20 km)

bbs_alpha <- list()

for (i in 1:length(radii)) {
	bbs_alpha[[i]] <- bbs_alphadiv %>%
		mutate(radius = radii[i]) %>%
		group_by(rteNo, radius) %>%
		summarize_at(c('richness', 'MPD_pa_z', 'MNTD_pa_z', 'MPDfunc_pa_z', 'MNTDfunc_pa_z'), funs(median(.[stop_bounds[i,1]:stop_bounds[i,2]][is.finite(.[stop_bounds[i,1]:stop_bounds[i,2]])]))) %>%
		ungroup 
}

bbs_alpha <- bbs_alpha %>% bind_rows %>% arrange(rteNo, radius)

#write.csv(bbs_alpha, file = '/mnt/research/nasabio/data/bbs/bbs_withinroute_alpha.csv', row.names = FALSE)	
write.csv(bbs_alpha, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/withinroute/bbs_withinroute_alpha_with1km.csv', row.names = FALSE)	

# Compile beta diversity: BBS within route
# ----------------------------------------

# It's already averaged by radius so we just need to put the slices together.

load('/mnt/research/nasabio/data/bbs/bbsworkspace_bystop_20072016.r')
library(dplyr)

radii <- c(1, 5, 10, 20) # in km

bbs_beta <- lapply(1:3089, function(i) {
	load(paste0('/mnt/research/nasabio/data/bbs/diversitywithinroute/beta_',i,'.r'))
	beta_div
})

#bbs_beta <- do.call(c, bbs_beta) %>% bind_rows 
bbs_beta <- bind_rows(bbs_beta)

route_ids <- unique(bbscov_oneyear$rteNo)

bbs_beta <- cbind(rteNo = rep(route_ids, each=length(radii)), bbs_beta)
bbs_beta <- bbs_beta[, apply(bbs_beta, 2, function(x) any(!is.na(x)))]

write.csv(bbs_beta, file = '/mnt/research/nasabio/data/bbs/biodiversity_CSVs/withinroute/bbs_withinroute_beta.csv', row.names = FALSE)	
