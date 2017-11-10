# Compile betapart for FIA PNW
# Using new coordinates
# 10 Nov 2017

fia_betapart <- list()
radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300) * 1000

# Only have up to 100 km right now.
radii <- radii[1:6]
nslices <- rep(250, 6)

for (r in 1:length(radii)) {
	print(r)
	list_r <- list()
	for (slice in 1:nslices[r]) {

		dat <- read.csv(paste0('/mnt/research/nasabio/data/fia/betaoutput/sliced2/betabaselga_',as.character(as.integer(radii[r])),'_',slice,'.csv'))
		list_r[[slice]] <- dat

	}

fia_betapart[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))

}

fia_betapart <- do.call('rbind', fia_betapart)

# Add identifying information
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords.csv')

gamma_div_cast <- cbind(plotmetadata, gamma_div_cast[, !names(gamma_div_cast) %in% 'plot'])
write.csv(gamma_div_cast, file = '/mnt/research/nasabio/data/fia/fia_gammadiv.csv')
