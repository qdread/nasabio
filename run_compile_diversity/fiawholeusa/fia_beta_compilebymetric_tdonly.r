# Make a sparse matrix for each of the metrics.
# Edited 30 March: need to store in a less memory intensive way. *was a full size matrix*
# Edited 28 Nov 2018: use SLURM IDs

# Function to take a vector and return a two column matrix
# First column is index, second column is value
col2sparse <- function(x) {
	notna <- !is.na(x)
	cbind(which(notna), x[notna])
}

column <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

n_plot <- 135032

metric_list <- list()
fp <- '/mnt/research/nasabio/data/fia/diversity/usa2018'

for (i in 1:n_plot) {
	load(file.path(fp, paste0('td100beta_', i, '.r')))
	metric_list[[i]] <- col2sparse(beta_div[,column])
	if(i%%1000==0) print(i)
	rm(beta_div) # Remove and do garbage collection.
	gc()
}

save(metric_list, file = file.path(fp, paste0('metric_betatd100_', column, '.r')))
