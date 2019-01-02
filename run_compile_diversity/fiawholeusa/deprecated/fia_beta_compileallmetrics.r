# Make a sparse matrix for ALL beta diversity metrics.

# Function to create sparse matrix from large matrix
# First column is the index, remaining columns are values
allmetrics2sparse <- function(x) {
	notna <- apply(x, 1, function(z) any(!is.na(z)))
	cbind(which(notna), x[notna, , drop = FALSE])
}

n_plot <- 119177

metric_list <- list()
fp <- '/mnt/research/nasabio/data/fia/diversity/usa2018'

for (i in 1:n_plot) {
	load(file.path(fp, paste0('beta_', i, '.r')))
	beta_div_i <- cbind(i, allmetrics2sparse(beta_div))
	beta_div_i <- rbind(beta_div_i, beta_div_i[, c(2, 1, 3:ncol(beta_div_i))])
	if (dim(beta_div_i)[2] == 23) metric_list[[i]] <- beta_div_i else metric_list[[i]] <- matrix(0, nrow=0, ncol=23)
	if(i%%1000==0) print(i)
	rm(beta_div) # Remove and do garbage collection.
	gc()
}

metric_list <- do.call(rbind, metric_list)
save(metric_list, file = file.path(fp, 'sparsematrixallbetametrics.r'))
