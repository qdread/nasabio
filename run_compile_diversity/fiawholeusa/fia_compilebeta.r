###################################
# Compile fia beta-diversity for whole USA

# Update 08 Sep: done again with new traits, now there are 10,000 slices.
# Update 03 Dec: done again for new FIA diversity (unfuzzed)
# Update 08 Jan 2018: do entire USA FIA beta-diversity (only taxonomic for now)


# Taxonomic only ----------------------------------------------------------


n_slices <- 250

fia_taxbetadiv <- list()

pb <- txtProgressBar(0, n_slices, style = 3)

for (i in 1:n_slices) {
	load(paste0('/mnt/research/nasabio/data/fia/diversity/usa/tdbeta_', i, '.r'))
  # Delete excess columns to save space.
  beta_div_list <- lapply(beta_div_list, '[', , 1:5)
	fia_taxbetadiv[[i]] <- beta_div_list
	setTxtProgressBar(pb, i)
}

close(pb)
fia_taxbetadiv <- do.call(c, fia_taxbetadiv)

save(fia_taxbetadiv, file = '/mnt/research/nasabio/data/fia/fia_betadivtd_listbypoint.r')
# Results in a list of pairwise matrices.

load('/mnt/research/nasabio/data/fia/fia_betadivtd_listbypoint.r')
fia_betadiv_list <- list()

for (divmetric in 1:ncol(fia_taxbetadiv[[1]])) {
  print(divmetric)
  list_i <- list()
  for (plot in 1:length(fia_betadiv)) {
    list_i[[plot]] <- fia_betadiv[[plot]][, divmetric]
  }
  fia_betadiv_list[[divmetric]] <- do.call('cbind', list_i)
}

names(fia_betadiv_list) <- dimnames(fia_betadiv[[1]])[[2]]

save(fia_betadiv_list, file = '/mnt/research/nasabio/data/fia/fia_betadivtdpdfd_listbymetric.r')


# All diversity types -----------------------------------------------------

n_slices <- 135174

fia_betadiv <- list()

pb <- txtProgressBar(0, n_slices, style = 3)

for (i in 1:n_slices) {
	load(paste0('/mnt/research/nasabio/data/fia/diversity/usa/beta_', i, '.r'))
	fia_betadiv[[i]] <- beta_div
	setTxtProgressBar(pb, i)
}

close(pb)

save(fia_betadiv, file = '/mnt/research/nasabio/data/fia/fia_usa_betadivtdpdfd_listbypoint.r')


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
