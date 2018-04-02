# Generate a single matrix 135174 x 135174 for each of the metrics.
# Edited 30 March: need to store in a less memory intensive way.

# Function to take a vector and return a two column matrix
# First column is index, second column is value
col2sparse <- function(x) {
	notna <- !is.na(x)
	cbind(which(notna), x[notna])
}

column <- as.numeric(Sys.getenv('PBS_ARRAYID'))

n_plot <- 135174

metric_list <- list()
fp <- '/mnt/research/nasabio/data/fia/diversity/usa'

for (i in 1:n_plot) {
	load(file.path(fp, paste0('beta_', i, '.r')))
	metric_list[[i]] <- col2sparse(beta_div[,column])
	if(i%%1000==0) print(i)
	rm(beta_div) # Remove and do garbage collection.
	gc()
}

#library(dplyr)
#metric_list <- do.call(rbind, metric_list)

save(metric_list, file = file.path(fp, paste0('metric_beta_', column, '.r')))
