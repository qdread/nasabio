###################################
# 23 Jun: compile FIA beta diversity lookup table

# Update 08 Sep: done again with new traits, now there are 10,000 slices.
# Update 03 Dec: done again for new FIA diversity (unfuzzed)

n_slices <- 10000

fia_betadiv <- list()

pb <- txtProgressBar(0, n_slices, style = 3)

for (i in 1:n_slices) {
	load(paste0('/mnt/research/nasabio/data/fia/diversity/unfuzzed/beta_', i, '.r'))
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
	fia_betadiv_list[[divmetric]] <- do.call('cbind', list_i)
}

names(fia_betadiv_list) <- dimnames(fia_betadiv[[1]])[[2]]

save(fia_betadiv_list, file = '/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_listbymetric.r')
