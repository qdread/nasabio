bbs_betadiv <- list()

radii <- c(1000,2000,3000,4000,5000)

for (r in 1:length(radii)) {
	list_r <- list()
	for (slice in 1:50) {

dat <- read.csv(paste0('/mnt/research/nasabio/data/bbs/betaoutput/td_',radii[r],'_',slice,'.csv'))
list_r[[slice]] <- dat

}

bbs_betadiv[[r]] <- data.frame(radius = radii[r], do.call('rbind', list_r))

}

# Combine with the covariates.

load('/mnt/research/nasabio/data/bbs/bbsworkspace.r')

bbs_betadiv <- cbind(bbscov, do.call('rbind', bbs_betadiv))
write.csv(bbs_betadiv, file = '/mnt/research/nasabio/data/bbs/bbs_beta_td.csv', row.names = FALSE)
