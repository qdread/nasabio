bbs_betadiv <- list()

radii <- c(1000,2000,3000,4000,5000)

for (r in 1:length(radii)) {
	list_r <- list()
	for (slice in 1:50) {

dat <- read.csv(paste0('/mnt/research/nasabio/data/bbs/betaoutput/td/td_',radii[r],'_',slice,'.csv'))
list_r[[slice]] <- dat

}

bbs_betadiv[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))

}

# Combine with the covariates.

load('/mnt/research/nasabio/data/bbs/bbsworkspace.r')

bbs_betadiv <- cbind(bbscov, do.call('rbind', bbs_betadiv))
write.csv(bbs_betadiv, file = '/mnt/research/nasabio/data/bbs/bbs_beta_td.csv', row.names = FALSE)

###################################################

# Compile BBS route-aggregated beta-diversity.

bbs_betadiv_route <- list()

radii <- c(50000, 75000, 100000)

for (r in 1:length(radii)) {
	list_r <- list()
	for (slice in 1:50) {

dat <- read.csv(paste0('/mnt/research/nasabio/data/bbs/betaoutput/byroute/route_tdpdfd_',as.character(as.integer(radii[r])),'_',slice,'.csv'))
list_r[[slice]] <- dat

}

bbs_betadiv_route[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))

}

load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
bbs_betadiv_route <- cbind(bbscov, do.call('rbind', bbs_betadiv_route))
write.csv(bbs_betadiv_route, file = '/mnt/research/nasabio/data/bbs/bbs_beta_byroute.csv', row.names = FALSE)


###############################
# 19 May: compile BBS alpha diversity.

bbs_alphadiv <- list()

for (i in 1:250) {
	load(paste0('/mnt/research/nasabio/data/bbs/diversity/alpha_',i,'.r'))
	bbs_alphadiv[[i]] <- alpha_div
}

bbs_alphadiv <- do.call('rbind', bbs_alphadiv)
write.csv(bbs_alphadiv, file = '/mnt/research/nasabio/data/bbs/bbs_alphadiv.csv', row.names = FALSE)