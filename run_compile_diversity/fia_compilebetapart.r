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

		dat <- read.csv(paste0('/mnt/research/nasabio/data/fia/diversity/unfuzzed/betapartfinal_',as.character(as.integer(radii[r])),'_',slice,'.csv'))
		list_r[[slice]] <- dat

	}

fia_betapart[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))

}



# Add identifying information
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords.csv')

fia_betapart <- cbind(plotmetadata[rep(1:nrow(plotmetadata), each = 9), ], do.call('rbind', fia_betapart))

write.csv(fia_betapart, file = '/mnt/research/nasabio/data/fia/fia_betapartfinal_to100.csv', row.names = FALSE)

# Edit: make this into a wide format data frame, similar to the alpha and gamma ones.
library(reshape2)
fia_betapart_wide <- dcast(fia_betapart, PLT_CN + STATECD + COUNTYCD + PLOT + radius + nneighb ~ index + diversity + partition)

write.csv(fia_betapart_wide, file = '/mnt/research/nasabio/data/fia/fia_betapartfinal_to100_wide.csv', row.names = FALSE)