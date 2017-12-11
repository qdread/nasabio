# Compile fia beta slices. Must submit job because it gets killed if run on dev node.

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
