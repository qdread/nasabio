# Compile BBS route-aggregated beta-diversity partitioning metrics (Baselga and Podani, taxonomic only)

bbs_betapart_route <- list()

radii <- c(50, 75, 100, 150, 200, 300) * 1000

for (r in 1:length(radii)) {
	list_r <- list()
	for (slice in 1:20) {

dat <- read.csv(paste0('/mnt/research/nasabio/data/bbs/betaoutput/byroute/route_betabaselga_',as.character(as.integer(radii[r])),'_',slice,'.csv'))
list_r[[slice]] <- dat

}

bbs_betapart_route[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))

}

load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
bbs_betapart_route <- cbind(bbscov[rep(1:nrow(bbscov), each = 20), ], do.call('rbind', bbs_betapart_route))
write.csv(bbs_betapart_route, file = '/mnt/research/nasabio/data/bbs/bbs_betapart_byroute.csv', row.names = FALSE)

