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

###############################
# 19 May: compile FIA alpha + gamma diversity.

fia_alphadiv <- list()

for (i in 1:250) {
	load(paste0('/mnt/research/nasabio/data/fia/diversity/alpha_',i,'.r'))
	fia_alphadiv[[i]] <- alpha_div
}

fia_alphadiv <- do.call('rbind', fia_alphadiv)
write.csv(fia_alphadiv, file = '/mnt/research/nasabio/data/fia/fia_alphadiv.csv', row.names = FALSE)

fia_gammadiv <- list()
rowidx <- round(seq(0,22532,length.out=251))

library(reshape2)

for (i in 1:250) {
	print(i)
	rowidxmin <- rowidx[i]+1
	rowidxmax <- rowidx[i+1]

	load(paste0('/mnt/research/nasabio/data/fia/diversity/gamma_',i,'.r'))
	gamma_div_melt <- melt(gamma_div[rowidxmin:rowidxmax,,], varnames = c('plot','radius','diversity'))
	gamma_div_cast <- dcast(gamma_div_melt, plot + radius ~ diversity)
	gamma_div_cast$radius <- as.numeric(substr(gamma_div_cast$radius, 3, nchar(as.character(gamma_div_cast$radius))))
	fia_gammadiv[[i]] <- gamma_div_cast
}

fia_gammadiv <- do.call('rbind', fia_gammadiv)

# Join with identifiers
load('/mnt/research/nasabio/data/fia/fiaworkspace2.r')
fia_gammadiv <- cbind(fiacoords[rep(1:nrow(fiacoords), each=14),], fia_gammadiv[,-1])


write.csv(fia_gammadiv, file = '/mnt/research/nasabio/data/fia/fia_gammadiv.csv', row.names = FALSE)

###################################
# 23 Jun: compile FIA beta diversity lookup table

# Update 08 Sep: done again with new traits, now there are 10,000 slices.

n_slices <- 10000

fia_betadiv <- list()

pb <- txtProgressBar(0, n_slices, style = 3)

for (i in 1:n_slices) {
	load(paste0('/mnt/research/nasabio/data/fia/diversity/beta_', i, '.r'))
	fia_betadiv[[i]] <- beta_div
	setTxtProgressBar(pb, i)
}

close(pb)
fia_betadiv <- do.call(c, fia_betadiv)

save(fia_betadiv, file = '/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_listbypoint.r')

# Results in a list of pairwise matrices.

# Go through and convert the list into a list of matrices that are essentially pairwise distance matrices.
# An array can be created for each year. route x route x beta-diversity metric. (or a list of equal size matrices)
# Then put those arrays into a list.
### 08 Sep update: Not done because it allocates too much RAM. Instead, let's get the needed values out on a case by case basis.

fia_betadiv_array <- array(NA, dim = c(nrow(fia_betadiv[[1]]), nrow(fia_betadiv[[1]]), ncol(fia_betadiv[[1]])))

for (divmetric in 1:ncol(fia_betadiv[[1]]))
	fia_betadiv_array[,,divmetric] <- do.call(cbind, lapply(fia_betadiv, '[', , divmetric))

fia_betadiv_array <- lapply(fia_betadiv_array, function(x) {dimnames(x)[[3]] <- dimnames(fia_betadiv[[1]])[[2]]; x})

save(fia_betadiv_array, file = '/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_array.r')

####
# Instead of the array, crate pairwise matrices in a list, probably will allocate less memory.
load('/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_listbypoint.r')
fia_betadiv_list <- list()

for (divmetric in 1:ncol(fia_betadiv[[1]])) {
	print(divmetric)
	list_i <- list()
	for (plot in 1:length(fia_betadiv)) {
		list_i[[plot]] <- fia_betadiv[[plot]][, divmetric]
	}
	fia_betadiv_list[[i]] <- do.call('cbind', list_i)
}

names(fia_betadiv_list) <- dimnames(fia_betadiv[[1]])[[2]]

save(fia_betadiv_list, file = '/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_listbymetric.r')

