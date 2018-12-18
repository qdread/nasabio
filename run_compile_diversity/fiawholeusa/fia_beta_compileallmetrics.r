# Make a sparse matrix for ALL beta diversity metrics.

# Function to take a vector and return a two column matrix
# First column is index, second column is value
allmetrics2sparse <- function(x) {
	notna <- apply(x, 1, function(z) any(!is.na(z)))
	cbind(which(notna), x[notna, ])
}

n_plot <- 119177

metric_list <- list()
fp <- '/mnt/research/nasabio/data/fia/diversity/usa2018'

for (i in 1:n_plot) {
	load(file.path(fp, paste0('beta_', i, '.r')))
	metric_list[[i]] <- cbind(i, allmetrics2sparse(beta_div))
	if(i%%1000==0) print(i)
	rm(beta_div) # Remove and do garbage collection.
	gc()
}

metric_list <- do.call(rbind, metric_list)
save(metric_list, file = file.path(fp, 'sparsematrixallbetametrics.r'))
