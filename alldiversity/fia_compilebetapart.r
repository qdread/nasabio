# Compile FIA beta-diversity partitioning metrics (Baselga and Podani, taxonomic only)

fia_betapart <- list()

radii <- c(1, 5, 7.5, 10, 20, 50, 75, 100, 150, 200, 300) * 1000
nslices <- c(rep(20, 9), 250, 250)

for (r in 1:length(radii)) {
	list_r <- list()
	for (slice in 1:nslices[r]) {

dat <- read.csv(paste0('/mnt/research/nasabio/data/fia/betaoutput/sliced/betabaselga_',as.character(as.integer(radii[r])),'_',slice,'.csv'))
list_r[[slice]] <- dat

}

fia_betapart[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))

}

load('/mnt/research/nasabio/data/fia/fiaworkspace2.r')
fia_betapart <- cbind(fiacoords[rep(1:nrow(fiacoords), each = 20), ], do.call('rbind', fia_betapart))
write.csv(fia_betapart, file = '/mnt/research/nasabio/data/fia/fia_betapart.csv', row.names = FALSE)

