# Compile the beta diversity from fia

fp <- '/mnt/research/nasabio'
load(file.path(fp, 'data/fia/fia_taxonomicbetadiv3.RData')
fia_betadiv3k <- data.frame(beta_td_shannon = fia_shannonbetadiv, beta_td_pairwise = fia_meanpairwisedissim, beta_pd_nearesttaxon_z = fia_phypairwise_z, beta_pd_pairwise_z = fia_phynt_z)
save(fia_betadiv3k, file = file.path(fp, 'data/fia/fia_betadiv3k.r'))

#########################

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/'
load(file.path(fp,'fia_diversitymetrics.RData'))
load(file.path(fp,'fia_betadiv3k.r'))

fia_betadiv3k <- cbind(fiaalbers@data, fia_betadiv3k)
plot_diversity <- left_join(plot_diversity, fia_betadiv3k)
names(plot_diversity)[13:16] <- paste0(names(plot_diversity)[13:16], '3000')

#########################

# Compile beta-diversity from longer runs.

load('/mnt/research/nasabio/data/fia/fiaworkspace.r') 

radii <- c(5000, 7500, 10000)

fia_betadiv <- list()

for (r in 1:length(radii)) {
	dat <- read.csv(paste0('/mnt/research/nasabio/data/fia/betaoutput/tdpdfd_',as.character(as.integer(radii[r])),'.csv'))
	fia_betadiv[[r]] <- data.frame(radius = radii[r], dat)
}

fia_betadiv <- do.call('rbind', fia_betadiv)

##########################

# Compile taxonomic-only beta-diversity from 100km run.

load('/mnt/research/nasabio/data/fia/fiaworkspace.r') 

fia_betadiv <- list()

for (i in 1:50) {
	dat <- read.csv(paste0('/mnt/research/nasabio/data/fia/betaoutput/sliced/td_100000_',i,'.csv'))
	fia_betadiv[[i]] <- dat
}

fia_betadiv <- do.call('rbind', fia_betadiv)

# add to the existing data frame.
fia_betatd <- read.csv('/mnt/research/nasabio/data/fia/fia_betatd.csv', stringsAsFactors=FALSE)

fia_betadiv <- cbind(fia_betatd[1:nrow(fia_betadiv), 1:6], radius = 1e5, fia_betadiv)
fia_betatd <- rbind(fia_betatd, fia_betadiv)
write.csv(fia_betatd, file = '/mnt/research/nasabio/data/fia/fia_betatd.csv', row.names = FALSE)