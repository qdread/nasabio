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


################################
# 22 May: compile taxonomic-only route-aggregated BBS beta diversity, by year.
# These are pairwise comparisons. It will generate a lookup table that I will go through and get median pairwise for each radius.


years <- 1997:2015
n_slices <- 13

bbs_betadivtd <- replicate(length(years), list()) # One list for each year.

combos <- expand.grid(1:n_slices, years)

for (i in 1:nrow(combos)) {
	load(paste0('/mnt/research/nasabio/data/bbs/diversity/tdbeta_', combos[i,2], '_', combos[i,1], '.r'))
	bbs_betadivtd[[which(years==combos[i,2])]] <- c(bbs_betadivtd[[which(years==combos[i,2])]], beta_div)
}

# Results in a list of lists. 19 years long, each one is a pairwise.

# Go through and convert the list into a list of matrices that are essentially pairwise distance matrices.
# An array can be created for each year. route x route x beta-diversity metric. (or a list of equal size matrices)
# Then put those arrays into a list.

bbs_betadivtd_arraylist <- lapply(bbs_betadivtd, function(x) array(NA, dim = c(nrow(x[[1]]), nrow(x[[1]]), ncol(x[[1]]))))

for (year in 1:length(bbs_betadivtd)) {
	for (divmetric in 1:ncol(bbs_betadivtd[[1]][[1]]))
	bbs_betadivtd_arraylist[[year]][,,divmetric] <- do.call(cbind, lapply(bbs_betadivtd[[year]], '[', , divmetric))

}

bbs_betadivtd_arraylist <- lapply(bbs_betadivtd_arraylist, function(x) {dimnames(x)[[3]] <- dimnames(bbs_betadivtd[[1]][[1]])[[2]]; x})
names(bbs_betadivtd_arraylist) <- paste0('year',years)

save(bbs_betadivtd_arraylist, file = '/mnt/research/nasabio/data/bbs/bbs_betadiv_arraylist.r')

################################
# 23 May: compile BBS gamma-diversity, by year
# Updated 01 June with latest year

years <- 1997:2016
n_slices <- 12
combos <- expand.grid(1:n_slices, years)

bbs_gammadiv <- replicate(length(years), list()) # One list for each year.

library(reshape2)

load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
yeartable <- table(bbscov$year)



# Combine arrays since I stupidly put wrong indices.
for (i in 1:length(years)) {
	print(i)
	rowidx <- round(seq(0,yeartable[which(names(yeartable)==years[i])],length.out=n_slices + 1))
	
	
	for (j in 1:n_slices) {
		rowidxmin <- rowidx[j]+1
		rowidxmax <- rowidx[j+1]
		load(paste0('/mnt/research/nasabio/data/bbs/diversity/gamma_', years[i], '_', j, '.r'))
		gamma_div_melt <- melt(gamma_div[rowidxmin:rowidxmax,,], varnames = c('rteNo','radius','diversity'))
		gamma_div_cast <- dcast(gamma_div_melt, rteNo + radius ~ diversity)
		gamma_div_cast$radius <- as.numeric(substr(gamma_div_cast$radius, 3, nchar(as.character(gamma_div_cast$radius))))
		bbs_gammadiv[[i]][[length(bbs_gammadiv[[i]]) + 1]] <- gamma_div_cast
		
	}
	
	bbs_gammadiv[[i]] <- do.call(rbind, bbs_gammadiv[[i]])
	cov_i <- subset(bbscov, year == years[i])
	bbs_gammadiv[[i]] <- cbind(cov_i[rep(1:nrow(cov_i), each=9),], bbs_gammadiv[[i]][,-1])
	
}

bbs_gammadiv <- do.call(rbind, bbs_gammadiv)
write.csv(bbs_gammadiv, file = '/mnt/research/nasabio/data/bbs/bbs_gammadiv.csv', row.names = FALSE)


